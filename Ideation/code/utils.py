import numpy as np
import torch
import dgl

def pbar(p=0, msg="", bar_len=20):
    msg = msg.ljust(50)
    block = int(round(bar_len * p))
    text = "\rProgress: [{}] {}% {}".format("\x1b[32m" + "="*(block-1) + ">" + "\033[0m" + "-" * (bar_len - block), round(p * 100, 2), msg)
    print(text, end="")
    if p == 1:
        print()

class AvgMeter:
    def __init__(self):
        self.reset()

    def reset(self):
        self.metrics = {}

    def add(self, batch_metrics):
        if self.metrics == {}:
            for key, value in batch_metrics.items():
                self.metrics[key] = [value]
        else:
            for key, value in batch_metrics.items():
                self.metrics[key].append(value)

    def get(self):
        return {key: np.mean(value) for key, value in self.metrics.items()}

    def msg(self):
        avg_metrics = {key: np.mean(value) for key, value in self.metrics.items()}
        return "".join(["[{}] {:.5f} ".format(key, value) for key, value in avg_metrics.items()])
    
def analyse(graph, actual_graph, vae_graph, edges, metapath):

    canonical_etypes = graph.canonical_etypes

    unique = []
    
    for i in canonical_etypes:
        src, r, dst = i
        unique.append(r)

    unique = sorted(np.unique(unique))

    edge_dict = {key:i for i, key in enumerate(unique)}

    nodes = {n_type : graph.num_nodes(n_type) for n_type in sorted(graph.ntypes)}

    num_nodes = graph.num_nodes()

    vae_canonical_etypes = vae_graph.canonical_etypes

    vae_unique = []
    
    for i in vae_canonical_etypes:
        src, r, dst = i
        vae_unique.append(r)

    vae_unique = sorted(np.unique(vae_unique))

    vae_nodes = {n_type : vae_graph.num_nodes(n_type) for n_type in sorted(vae_graph.ntypes)}

    vae_unique = sorted(np.unique(vae_unique))

    vae_edge_dict = {key:i for i, key in enumerate(vae_unique)}

    relevant = [edge_dict[i] for i in vae_unique]

    removed = []

    for ind, etype in enumerate(metapath):

        temp = (torch.stack(actual_graph.edges(etype = etype), 0)).T[edges[etype]]
        b = torch.ones((temp.shape[0], 1))*ind
        temp = torch.cat([b, temp], -1)
        removed.append(temp)

    removed = torch.cat(removed, 0).long()    

    link_labels = torch.cat(([vae_graph.adj(etype = i).to_dense().view(-1) for i in metapath]))

    new = dgl.edge_type_subgraph(actual_graph, metapath)

    actual_link_labels = torch.cat(([new.adj(etype = i).to_dense().view(-1) for i in new.canonical_etypes]))
        
    return {"nodes" : nodes, "edges" : unique, "edge_dict" : edge_dict, 
            "canonical_etypes" : canonical_etypes, "num_nodes" : num_nodes, 
            "num_edges" : len(unique), "removed" : removed, 
            "link_labels" : link_labels, 'relevant' : relevant, "vae_nodes" : vae_nodes, "vae_edges" : vae_unique, 
            "vae_canonical_etypes" : vae_canonical_etypes, "vae_edge_dict" : vae_edge_dict, "actual_link_labels":actual_link_labels}
        


