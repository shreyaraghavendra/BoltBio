import dgl
import torch
import dgl.function as fn
import torch.nn as nn


def u_mul_e(emb):

    def edge_message(edges):

        feat =  edges.src['feat']

        return {'m' : feat * emb.repeat(feat.shape[0], 1)}
    
    return edge_message

def u_mul_e_sub_v(emb):

    def edge_message(edges):

        feat =  edges.src['feat']

        return {'m' : edges.dst['feat'] - feat * emb.repeat(feat.shape[0], 1)}
    
    return edge_message

class GConvAttn(nn.Module):

    def __init__(self, dim, hid_dim, edge_dict):

        super(GConvAttn, self).__init__()

        self.fc_q = nn.Linear(dim, hid_dim, bias = False)
        self.fc_k = nn.Linear(dim, hid_dim, bias = False)
        self.fc_v = nn.Linear(dim, hid_dim, bias = False)
        self.attn_fc = nn.Linear(hid_dim, hid_dim) 
        self.post_fc = nn.Linear(hid_dim, dim)

        self.edge_dict = edge_dict

    def forward(self, g, emb):

        with g.local_scope():

            funcs = {}

            for c_etype in g.canonical_etypes:
                
                src, etype, dst = c_etype
                funcs[c_etype] = (u_mul_e(emb[self.edge_dict[etype]]), fn.mean('m', 'key'))

            g.multi_update_all(funcs, "stack")

            for n_type in g.ntypes:
                if 'key' in g.nodes[n_type].data.keys():
                    k = g.nodes[n_type].data['key']
                    del g.nodes[n_type].data['key']

                    q = g.nodes[n_type].data['feat']

                    res = q
                    
                    q = self.fc_q(q)[:,None,:].repeat( 1, k.shape[1], 1)
                    k, v = self.fc_k(k), self.fc_v(k)

                    a = torch.softmax(self.attn_fc(q - k), dim=1)
                    o = (a*(v)).sum(dim = 1)
                    o = self.post_fc(o) + res

                    g.nodes[n_type].data['feat'] = o

            return g.ndata['feat']

class GCN_node(nn.Module):

    def __init__(self, in_dim, out_dim, edge_dict):

        super(GCN_node, self).__init__()

        self.mlp = nn.Sequential(nn.Linear(in_dim, out_dim, bias = False), nn.BatchNorm1d(out_dim), nn.ReLU(), nn.Linear(out_dim, out_dim))

        self.edge_dict = edge_dict
    
    def forward(self, g, emb):

        with g.local_scope():

            funcs = {}

            for c_etype in g.canonical_etypes:
                
                src, etype, dst = c_etype

                funcs[c_etype] = (u_mul_e_sub_v(emb[self.edge_dict[etype]]), fn.mean('m', 'k'))

            g.multi_update_all(funcs, "sum")

            keys, k = list(g.ndata['k'].keys()), list(g.ndata['k'].values())
            lens = [i.shape[0] for i in k]
            k = self.mlp(torch.cat(k, 0))
            feat = torch.cat([g.nodes[i].data['feat'] for i in keys], 0)
            feat = feat + k
            feat = feat.split(lens, 0)
            del k

            for key, f in zip(keys, feat):

                g.nodes[key].data['feat'] = f

            return g.ndata['feat']

class MLP(nn.Module):

    def __init__(self, in_dim, out_dim):

        super(MLP, self).__init__()

        self.mlp_nodes = nn.Sequential(nn.Linear(in_dim, out_dim, bias = False), nn.BatchNorm1d(out_dim), nn.ReLU(), nn.Linear(out_dim, out_dim))

        self.mlp_e = nn.Sequential(nn.Linear(in_dim, out_dim, bias = False), nn.BatchNorm1d(out_dim), nn.ReLU(), nn.Linear(out_dim, out_dim))

    def forward(self, g, emb):

        with g.local_scope():

            keys, feat = list(g.ndata['feat'].keys()), list(g.ndata['feat'].values())
            lens = [i.shape[0] for i in feat]
            feat = self.mlp_nodes(torch.cat(feat, 0)).split(lens, 0)

            # for i, f in zip(keys, feat):
              
            #     g.nodes[i].data['feat'] = f

            g.ndata['feat'] = {key:f for key, f in zip(keys, feat)}                

            emb = self.mlp_e(emb)

            return g.ndata['feat'], emb

class Encoder(nn.Module):

    def __init__(self, num_nodes, num_edges, in_dim, hid_dim, edge_dict, **kwargs):

        super(Encoder, self).__init__()

        self.n_emb = nn.Embedding(num_nodes, in_dim)

        self.e_emb = nn.Parameter(torch.randn((num_edges, in_dim)).float())

        self.num_nodes = num_nodes
        
        self.num_edges = num_edges

        self.mlp1 = GCN_node(in_dim = in_dim, out_dim = in_dim, edge_dict = edge_dict)

        self.attn1 = GConvAttn(in_dim, hid_dim, edge_dict)

        self.mlp2 = GCN_node(in_dim = in_dim, out_dim = in_dim, edge_dict = edge_dict)

        self.attn2 = GConvAttn(in_dim, hid_dim, edge_dict)

    def forward(self, g):

        with g.local_scope():

            g.ndata['feat'] = {ntype : self.n_emb(g.nodes[ntype].data['ind']).float() for ntype in g.ntypes}

            g.ndata['feat'] = self.mlp1(g, self.e_emb)

            g.ndata['feat'] = self.attn1(g, self.e_emb)

            g.ndata['feat'] = self.mlp2(g, self.e_emb)

            g.ndata['feat'] = self.attn2(g, self.e_emb)
 
            return g.ndata['feat'], self.e_emb

class parameterize(nn.Module):

    def __init__(self,):

        super(parameterize, self).__init__()

    def forward(self, g, emb):

        with g.local_scope():

            keys, feat = list(g.ndata['feat'].keys()), list(g.ndata['feat'].values())
            lens = [i.shape[0] for i in feat]
            mean, logvar = torch.cat(feat, 0).chunk(2, -1)
            feat = torch.normal(mean, torch.exp(logvar) * 0.5).split(lens, 0)

            mean, logvar = emb.chunk(2, -1)
            emb = torch.normal(mean, torch.exp(logvar) * 0.5)
            
            g.ndata['feat'] = {i: f for i, f in zip(keys, feat)}

            return g.ndata['feat'], emb

class VAE(nn.Module):

    def __init__(self, in_dim, dims, encoder, edge_dict, **kwargs):
        
        super(VAE, self).__init__()

        self.mlps = nn.ModuleList()
        self.attn = nn.ModuleList()

        self.encoder = encoder

        inp_dim = in_dim

        for i in dims:

            mlp = MLP(inp_dim, i)
            attn = GConvAttn(i, i*2, edge_dict)
            self.mlps.append(mlp)
            self.attn.append(attn)
            inp_dim = i

        self.enc_mlp = MLP(inp_dim, inp_dim * 2)

        self.dec_mlps = nn.ModuleList()
        self.dec_attn = nn.ModuleList()

        self.parameterize = parameterize()
        
        dims.insert(0, in_dim)

        for i in dims[::-1][1:]:

            mlp = MLP(inp_dim * 2, i)
            attn = GConvAttn(i, i*2, edge_dict)
            self.dec_mlps.append(mlp)
            self.dec_attn.append(attn)
            inp_dim = i
    
    def forward(self, g):

        with g.local_scope():
            g.ndata['feat'], e_emb = self.encoder(g)

            

            encoder = [g.ndata['feat'], e_emb]
            
            feats = []
            embs = []

            for ind, (layer1, layer2) in enumerate(zip(self.mlps, self.attn)):
                g.ndata['feat'], e_emb = layer1(g, e_emb)
                feat = layer2(g, e_emb)
                feats.append(feat)
                embs.append(e_emb)
                g.ndata['feat'] = feat
            
            g.ndata['feat'], e_emb = self.enc_mlp(g, e_emb) 
            
            g.ndata['feat'], e_emb = self.parameterize(g, e_emb)
            
            save = (torch.cat([g.ndata['feat'][i] for i in g.ndata['feat'].keys()], 0), e_emb)
            

            for ind, (layer1, layer2) in enumerate(zip(self.dec_mlps, self.dec_attn)):
                g.ndata['feat'] = {key : torch.cat([g.nodes[key].data['feat'], feats[-ind - 1][key]], -1) for key in g.ndata['feat'].keys() }
                e_emb = torch.cat((e_emb, embs[-ind - 1]), -1)
                g.ndata['feat'], e_emb = layer1(g, e_emb)
                g.ndata['feat'] = layer2(g, e_emb)
            

            return encoder, save, g.ndata['feat'], e_emb

if __name__ == '__main__':
    graph = {("A", "1", "C") : (torch.tensor([0, 1]), torch.tensor([1, 2])), 
             ("C", "2", "B") : (torch.tensor([0, 1]), torch.tensor([1, 3])), 
             ("A", "2", "C") : (torch.tensor([0, 1]), torch.tensor([1, 2,]))}

    graph = dgl.heterograph(graph)

    for etype in graph.canonical_etypes:
        graph.edges[etype].data['ind'] = torch.randint(0, graph.num_edges(), (graph.num_edges(etype), ))
    for ntype in graph.ntypes:
        graph.nodes[ntype].data['ind'] = torch.randint(0, graph.num_nodes(), (graph.num_nodes(ntype), ))

    m = Encoder(graph.num_nodes(), graph.num_edges(), 5, 3)

    m(graph)

    for i in range(100):
        m(graph)




