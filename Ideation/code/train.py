import argparse
import yaml
import datetime
import os
import torch
import torch.nn.functional as F
import wandb
import random
import numpy as np
import utils
import graph
import dgl
import models
from dgl import DGLError
from utils import pbar
import torchmetrics
import pickle as pkl

class EmbCriterion():

    def __init__(self, nodes, edges, canonical_etypes, **kwargs):

        temp = []

        for ind, i in enumerate(nodes):
            temp.append(torch.ones((nodes[i], )) * ind) 

        binary = []

        for i in nodes:
            for r in edges:
                for j in nodes:
                    if (i, r, j) in canonical_etypes:
                        binary.append(1)
                    else:
                        binary.append(0)
        
        self.bce_labels = torch.tensor(binary).float()


        self.label = torch.cat(temp, 0).long()

        self.contrast_loss = torch.nn.CrossEntropyLoss()
        self.bce = torch.nn.BCEWithLogitsLoss()

    def __call__(self, n_feat, e_feat, g):

        with g.local_scope():

            temp = []
            temp2 = []
            lens = []
            for i in sorted(g.ntypes):
                feat = n_feat[i]
                temp.append(feat)
                lens.append(feat.shape[0])
                temp2.append(feat.mean(0, keepdims=True))

            means = torch.cat(temp2, 0)
            del temp2
            feats = torch.cat(temp, 0)
            del temp

            logits = feats @ means.permute(1, 0)
            label = self.label.to(logits.device)

            contrast_loss = self.contrast_loss(logits.float(), label)

            dim = feats.shape[-1]

            logits = (((means[:, None, :] * e_feat[None, :, :])).view(-1, dim) @ means.permute(1, 0) ).view(-1, )

            label = self.bce_labels.to(logits.device)

            bce_loss = self.bce(logits.float(), label)

            total_loss = bce_loss + contrast_loss

            return {'loss' : total_loss, 'node_loss' : contrast_loss, 'relation_loss' : bce_loss}

class LinkCriterion():

    def __init__(self, link_labels, relevant, vae_nodes, vae_edges, vae_canonical_etypes,  **kwargs):

        self.emb_crit = EmbCriterion(vae_nodes, vae_edges, vae_canonical_etypes)

        self.labels = link_labels.float()
        self.relevant = relevant


        self.mse = torch.nn.MSELoss()
        self.bce_loss = torch.nn.BCEWithLogitsLoss()
        self.kl = torch.nn.KLDivLoss(log_target=True)
        self.f1 = torchmetrics.F1(threshold=0.0505)
        self.roc = torchmetrics.AUROC()

    def __call__(self, gt, encoder, n_feat, e_feat, graph):

        e_feat = e_feat[self.relevant]
        
        loss_stats = self.emb_crit(n_feat, e_feat, graph) 

     

        n_types, feat = list(n_feat.keys()), list(n_feat.values())
        lens = [i.shape[0] for i in feat] 

        #feat = torch.cat(feat, 0)
        reconstruction_loss = self.mse(torch.cat(feat, 0), gt[0]) + self.mse(gt[1][self.relevant], e_feat)

        del gt
        

        src, dst = feat
        del feat
        dim = src.shape[-1]
       
        logits = (((src[:, None, :] * e_feat[None, :, :])).view(-1, dim) @ dst.permute(1, 0) ).view(lens[0], len(e_feat), lens[1],).permute(1, 0, 2).flatten()
        labels = self.labels.to(logits.device)

        link_loss = self.bce_loss(logits, labels)

        f1 = self.f1(logits.cpu(), labels.long().cpu())
        auroc = self.roc(logits.cpu(), labels.long().cpu())

        del logits, labels

        nodes, edges = encoder

        edges = edges[self.relevant]

        kl_loss = self.kl(nodes, torch.normal(0, 1, nodes.shape).to(nodes.device)) + self.kl(edges, torch.normal(0, 1, edges.shape).to(edges.device))

        loss = kl_loss + reconstruction_loss + 0.01 * loss_stats['loss'] + link_loss

        loss_stats['kl_loss'] = kl_loss
        loss_stats['reconstruction_loss'] = reconstruction_loss
        loss_stats['link_loss'] = link_loss
        loss_stats['loss'] = loss
        loss_stats['f1'] = f1
        loss_stats['auroc'] = auroc

        return loss_stats

class decode_output():
    
    def __init__(self, removed, vae_edge_dict, names, identifiers, actual_link_labels, thresh, relevant,  **kwargs):

        self.removed = removed
        self.edge_dict = vae_edge_dict
        self.names = names
        self.identifiers = identifiers
        self.thresh = thresh
        self.labels = actual_link_labels.long()
        self.accuracy = torchmetrics.Accuracy()
        self.f1 = torchmetrics.F1()
        self.relevant = relevant

    def __call__(self, n_feat, e_feat, graph):

        n_types, feat = list(n_feat.keys()), list(n_feat.values())
        lens = [i.shape[0] for i in feat] 

        e_feat = e_feat[self.relevant]

        src, dst = feat
        del feat

        dim = src.shape[-1]

        pred = (torch.sigmoid((((src[:, None, :] * e_feat[None, :, :])).view(-1, dim) @ dst.permute(1, 0) ).view(lens[0], len(e_feat), lens[1]).permute(1, 0, 2)) > self.thresh).long().cpu()

        ids = torch.stack(torch.where(pred == 1), 0).T

        pred = pred.flatten()

        accuracy = self.accuracy(pred, self.labels)
        f1 = self.f1(pred, self.labels)

        c = 0

        for i in self.removed:
            if i in ids:
                c += 1

        new_accuracy = c/len(self.removed)

        data = []

        for i in ids:
            r, n1, n2 = i
            r = int(r)
            n1 = int(n1)
            n2 = int(n2)
            r_name = list(self.edge_dict.keys())[r]
            node_p = 0
            pos = n1
            ind = int(graph.nodes[n_types[node_p]].data['ind'][pos])
            name1 = self.names[ind]
            identifiers1 = self.identifiers[ind]
            node_p = 1
            pos = n2
            ind = int(graph.nodes[n_types[node_p]].data['ind'][pos])
            name2 = self.names[ind]
            identifiers2 = self.identifiers[ind]
            temp = [r_name, name1, name2, identifiers1, identifiers2]
            data += [temp]

        return {'f1' : f1, 'accuracy' : accuracy, 'new_accuracy' : torch.tensor([new_accuracy]), 'data' : data}

class Trainer():

    def __init__(self, config, pretrain = True, reinit = False):

        self.config = config
        self.out_dir = config['output']
        self.device = torch.device('cuda' if torch.cuda.is_available else 'cpu')

        try:
            self.graph, self.actual_graph, self.sub_graph = dgl.load_graphs(f'{config["graph"]["save"]}/graph.bin')[0]

        except DGLError:

            graph.create_graph(**config['graph'])

            self.graph, self.actual_graph, self.sub_graph = dgl.load_graphs(f'{config["graph"]["save"]}/graph.bin')[0]

        with open(f'{config["graph"]["save"]}/metapaths.pickle', 'rb') as f:
            metapaths = pkl.load(f)
        with open(f'{config["graph"]["save"]}/names.pickle', 'rb') as f:
            names = pkl.load(f)
        with open(f'{config["graph"]["save"]}/identifiers.pickle', 'rb') as f:
            identifiers = pkl.load(f)
        with open(f'{config["graph"]["save"]}/edges.pickle', 'rb') as f:
            edges = pkl.load(f)
        graph_spec = utils.analyse(self.graph, self.actual_graph, self.sub_graph, edges, metapaths)
        
        
        model = models.Encoder(**self.config['encoder'], **graph_spec)

        if pretrain:
            self.model = model.to(self.device)
            self.embed_crit = EmbCriterion(**graph_spec)
        else:
            model.load_state_dict(torch.load(config['encoder_save']))
            self.model = models.VAE(**config['vae'], **graph_spec, encoder = model).to(self.device)
            self.link_crit = LinkCriterion(**graph_spec)
            self.decode = decode_output(**graph_spec, names = names, identifiers = identifiers, **config['decode'])
        
        self.main_thread = True
        
        self.optim = torch.optim.Adam(self.model.parameters(), lr=self.config['lr'], weight_decay=0, betas=(0.9, 0.999))
        self.sched = torch.optim.lr_scheduler.StepLR(self.optim, step_size=40, gamma=0.7)
    
        if self.main_thread:
            if pretrain:
                if self.config["wandb"]:
                    run = wandb.init(name = self.config["name"])
                    config['wandb_url'] = run.get_url()

            os.makedirs(config['output'], exist_ok=True)

            with open(os.path.join(config["output"], 'config.yaml'), 'w') as outfile:
                yaml.dump(self.config, outfile)

            print(f"Configuration")
            print("--------------------------------------------")
            print(yaml.dump(self.config))
            print("--------------------------------------------")
            print(f"Model parameters: {sum(p.numel() for p in self.model.parameters())/1e6}M")
        
        if os.path.exists(os.path.join(self.config["output"], "last.ckpt")) and not reinit:

            ckpt = torch.load(os.path.join(self.config["output"], "last.ckpt"))
        
            self.model.load_state_dict(ckpt["model"])
            self.optim.load_state_dict(ckpt["optim"])
            self.sched.load_state_dict(ckpt["sched"])
            self.start_epoch = ckpt["epoch"] + 1
            if self.main_thread:
                print(f"Loaded checkpoint, continuing from {self.start_epoch} epochs...")   
        else:
            self.start_epoch = 0
            self.logs = {"train": [], "val": []}
            if self.main_thread:
                print(f"Checkpoint not found, starting fresh...")

        self.train_steps = self.start_epoch 
        self.metric_meter = utils.AvgMeter()
        self.log_wandb = self.config['wandb']
        self.pretrain = pretrain

    def train_embed(self):
        self.model.train()
        self.metric_meter.reset()

        self.graph = self.graph.to(self.device)

        for i in range(self.config['pretrain_steps']):
        
            n_feat, e_feat = self.model(self.graph)

            loss_stats = self.embed_crit(n_feat, e_feat, self.graph)

            loss = loss_stats['loss']

            self.optim.zero_grad()
            loss.backward()
            self.optim.step()


            metrics = {f'train {s}' : v.item() for s,v  in loss_stats.items()}

            self.metric_meter.add(metrics)

            if self.main_thread:
                if self.log_wandb:
                    metrics["train step"] = self.train_steps
                    wandb.log(metrics)
            pbar(i/self.config['steps'], msg=self.metric_meter.msg())

            self.train_steps += 1

        pbar(1, msg=self.metric_meter.msg())
        self.sched.step()

    def train_link(self):

        self.model.train()
        self.metric_meter.reset()
        self.sub_graph = self.sub_graph.to(self.device)

        for i in range(self.config['steps']):
        
            encoder, latent, n_feat, e_feat = self.model(self.sub_graph)

            loss_stats = self.link_crit(encoder, latent, n_feat, e_feat, self.sub_graph)

            loss = loss_stats['loss']

            self.optim.zero_grad()
            loss.backward()
            self.optim.step()

            metrics = {f'train {s}' : v.item() for s,v  in loss_stats.items()}

            self.metric_meter.add(metrics)

            if self.main_thread:
                if self.log_wandb:
                    metrics["train step"] = self.train_steps
                    wandb.log(metrics)
            pbar(i/self.config['steps'], msg=self.metric_meter.msg())

            self.train_steps += 1

        pbar(1, msg=self.metric_meter.msg())
        self.sched.step()

    @torch.no_grad()
    def val_link(self):

        self.model.eval()
        self.metric_meter.reset()
        self.sub_graph = self.sub_graph.to(self.device)
        
        _, _, n_feat, e_feat = self.model(self.sub_graph)

        loss_stats = self.decode(n_feat, e_feat, self.graph)

        metrics = {f'val {s}' : v.item() for s,v  in loss_stats.items() if s != 'data'}

        self.metric_meter.add(metrics)

        if self.main_thread:
            if self.log_wandb:
                wandb.log(metrics)

            pbar(1, msg=self.metric_meter.msg())

        return loss_stats['data']
        
    def run(self):
        best_train_loss, best_val_accuracy = float("inf"), 0

        self.metric_meter.reset()

        n_epochs = self.config['pretrain_epochs'] if self.pretrain else self.config['epochs']
        for epoch in range(self.start_epoch, n_epochs):

            if self.main_thread:
                print(f"Epoch: {epoch}")
                print("---------------")

            if self.pretrain:
                self.train_embed()
            else:
                self.train_link()

            if self.main_thread:
                train_loss = self.metric_meter.get()["train loss"]

                if train_loss < best_train_loss:
                    print(
                        "\x1b[34m"
                        + f"train loss improved from {round(best_train_loss, 5)} to {round(train_loss, 5)}"
                        + "\033[0m"
                    )
                    best_train_loss = train_loss
                
                if self.pretrain:

                    torch.save(self.model.state_dict(),
                                    os.path.join(self.out_dir, "best_encoder.ckpt"),
                                )

                    if self.log_wandb:
                        wandb.log(
                            {
                                "best_loss": best_train_loss,
                            }
                        )
                else:
                    data = self.val_link()

                    val_accuracy = self.metric_meter.get()["val accuracy"]

                    if self.log_wandb:
                        wandb.log({"predictions" : wandb.Table(data = data, columns = ['edge_type', 'node_1_name', 'node_2_name', 'node_1_identifier', 'node_2_identifier'])})
                        
                        
                        wandb.log(
                                {
                                    "epoch": epoch,
                                    "train": train_loss,
                                    "lr": self.optim.param_groups[0]["lr"],
                                }
                            )

                    if val_accuracy > best_val_accuracy:
                        print(
                            "\x1b[33m"
                            + f"val accuracy improved from {round(best_val_accuracy, 5)} to {round(val_accuracy, 5)}"
                            + "\033[0m"
                        )
                        best_val_accuracy = val_accuracy
                        if self.log_wandb:
                            wandb.log(
                                {
                                    "best_accuracy": self.metric_meter.get()["val new_accuracy"],
                                }
                            )


                        torch.save(
                                self.model.state_dict(),
                                    os.path.join(self.out_dir, "best.ckpt"),
                                )
                torch.save(
                    {
                        "model": self.model.state_dict(),
                        "sched": self.sched.state_dict(),
                        "optim": self.optim.state_dict(),
                        "epoch": epoch,
                    },
                    os.path.join(self.out_dir, "last.ckpt"),
                )
            

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-c", "--config", 
                         type = str,
                         required = True,
                         help = 'path to config file')
    
    parser.add_argument("-o", "--output",
                        type = str, 
                        default = f"output/{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')}",
                        help = "path to output file"
                        )

    parser.add_argument("-w","--wandb", 
                        action = "store_true", 
                        help = 'wandb logging')

    parser.add_argument("-n", "--name", 
                        type = str,
                        default = f"{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')}", 
                        help = 'name of wandb run')

    args = parser.parse_args()
    
    try:
        config = yaml.safe_load(open(args.config, "r"))
    except:
        raise ValueError(f"Incorrect path to config file : {args.config}")

    if not config.get('output', False):
        config['output'] = args.output

    if args.wandb:
        if not config.get('name', False):
            config['name'] = args.name
    
    config['wandb'] = args.wandb

    trainer = Trainer(config, pretrain = True) 
    trainer.run()
    config['encoder_save'] = f"{config['output']}/best_encoder.ckpt"
    trainer = Trainer(config, pretrain = False, reinit = True) 
    trainer.run()
    
