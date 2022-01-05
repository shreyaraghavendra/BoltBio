import json
import dgl
from tqdm import tqdm
import torch
import numpy as np
import os
import pickle as pkl
from copy import deepcopy

def create_graph(path, save, metanodes, train_split, reverse = False):

    with open(path, 'r') as f:
        file = json.load(f)

    nodes = dict()
    names = dict()
    identifiers = dict()
    indx = dict()
    ind = 0
    for i in file['nodes']:

        _type = i['kind']
        node = i['identifier']
        name = i['name']
        

        if _type in nodes:
            type_nodes = nodes[_type]
            type_nodes.append(node)
            nodes[_type] = type_nodes
            type_indx = indx[_type]
            type_indx.append(ind)
            indx[_type] = type_indx
        else:
            nodes[_type] = [node]
            names[_type] = [name]
            indx[_type] = [ind]

        names[ind] = name
        identifiers[ind] = node

        ind += 1
            
    graph = dict()

    eind = {}

    ind = 0

    for edge in tqdm(file['edges']):

        node_type, nid = edge['source_id']
        target_type, tid = edge['target_id']

        
        
        nid = nodes[node_type].index(nid)
        tid = nodes[target_type].index(tid)

        relation = edge['kind']

        edge_type = (node_type, relation, target_type)

        if edge_type in graph:
            _nodes, _targets = graph[edge_type] 
            _nodes.append(nid)
            _targets.append(tid)
            graph[edge_type] = (_nodes, _targets)
            inds = eind[edge_type]
            inds.append(ind)
            eind[edge_type] = inds
        else:
            graph[edge_type] = ([nid], [tid])
            eind[edge_type] = [ind]

        if nid != tid and reverse:

            nid, tid = tid, nid
            node_type, target_type = target_type, node_type

            edge_type = (node_type, relation, target_type)

            if edge_type in graph:
                _nodes, _targets = graph[edge_type] 
                _nodes.append(nid)
                _targets.append(tid)
                graph[edge_type] = (_nodes, _targets)
                inds = eind[edge_type]
                inds.append(ind)
                eind[edge_type] = inds
            else:
                graph[edge_type] = ([nid], [tid])
                eind[edge_type] = [ind]

        ind += 1
        
    eind = {key : torch.tensor(value) for key, value in eind.items()} 

    graph = dgl.heterograph(graph)

    for _type in graph.ntypes:

        graph.nodes[_type].data['ind'] = torch.tensor(indx[_type])

    og_graph = deepcopy(graph)

    metapaths = []

    for etype in graph.canonical_etypes:
        src, r, dst = etype
        p = True

        for edge in metanodes:
            for node_types in edge:
                if node_types not in [src, dst]:
                    p = False
                    break  
        if p:
            metapaths.append(etype)

        
    nums = [graph.num_edges(p) for p in metapaths]
    p = [num/sum(nums) for num in nums]
    train_size = int(sum(nums) * train_split)
    val_size = sum(nums) - train_size

    dist = np.unique(np.random.choice(len(metapaths),(val_size,), True, p = p ), return_counts = True)[-1]

    sub_edge = {p : np.random.choice(n, (s,), replace = False) for p, n, s in zip(metapaths, nums, dist)}

    for p in metapaths:
        graph.remove_edges(eids = sub_edge[p], etype = p)

    vae_graph = dgl.edge_type_subgraph(graph, etypes = metapaths)

    os.makedirs(f'{save}', exist_ok=True)

    dgl.save_graphs(f'{save}/graph.bin', [graph, og_graph, vae_graph])


    with open(f'{save}/metapaths.pickle', 'wb') as f:
        pkl.dump(metapaths, f)
    with open(f'{save}/identifiers.pickle', 'wb') as f:
        pkl.dump(identifiers, f)
    with open(f'{save}/names.pickle', 'wb') as f:
        pkl.dump(names, f)
    with open(f'{save}/edges.pickle', 'wb') as f:
        pkl.dump(sub_edge, f)
