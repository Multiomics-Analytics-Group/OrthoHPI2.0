import os
import utils
import pandas as pd
import networkx as nx
import json
from joblib import Parallel, delayed


def get_clusters_from_config(config_file):
    """"clusters": [{ "key": string, "color": string, "clusterLabel": string }], #key:tax id, color:hex, clusterLabel: name"""
    clusters = []
    hosts = utils.read_config(filepath=config_file, field='hosts')
    parasites = utils.read_config(filepath=config_file, field='parasites')
    for parasite in parasites:
        p = parasites[parasite]
        clusters.append({"key": parasite, "color":p['color'], "clusterLabel":p['label']})
    
    for host in hosts:
        h = hosts[host]
        clusters.append({"key": host, "color":h['color'], "clusterLabel":h['label']})

    return clusters

def generate_cytoscape_network(edges_file_path, proteins, tissues, config_file, output_dir_path, n=2):
    """
    Generates a file with the full network in cytoscape json format:
    {"nodes": [{"key": string, #identifier
                "label": string, #name shown
                "tag": string, #tissue id
                "URL": string, #linkout
                "cluster": string, #tax id
                "x": float, # position x
                "y": float, # position y
                "score": float}, # node size (i.e centrality)],
    "edges": [[string, string, float]], #source key, target key, score
    "clusters": [{ "key": string, "color": string, "clusterLabel": string }], #key:tax id, color:hex, clusterLabel: name
    "tags":[{ "key": string, "image": string }]} #key:tissue id image:path
    
    :param str edges_file_path: path to file with edges and their attributes
    :param dict proteins: dictionary of valid proteins (key: protein id, value: protein name)
    :param dict tissues: dictionary of relevant tissues (key: BTO, value: tissue name)
    :param str config_file: path to config_file
    :param str output_dir: path to output directory where to store the graphs
    :param int n: number of jobs to generate host-parasite graphs in parallel
    """
    clusters = get_clusters_from_config(config_file)
    tags = [{"key": "Parasite protein", "image": "concept.svg"}, {"key":"Human protein", "image": "person.svg" }]
    
    with open(edges_file_path, 'r') as f:
        data = pd.read_csv(f, sep='\t')
        
    parasites = utils.read_config(filepath=config_file, field='parasites')
    hosts = utils.read_config(filepath=config_file, field='hosts')
    Parallel(n_jobs=n)(delayed(build_graph)(data, parasite, hosts, proteins, tissues, clusters, tags, output_dir_path) for parasite in parasites)
    
    
def generate_common_cytoscape_network(edges_file_path, proteins, tissues, config_file, output_dir_path, min_common=2):
    clusters = get_clusters_from_config(config_file)
    tags = [{"key": "Parasite protein", "image": "concept.svg"}, {"key":"Human protein", "image": "person.svg" }]
    with open(edges_file_path, 'r') as f:
        data = pd.read_csv(f, sep='\t')
        
    hosts = utils.read_config(filepath=config_file, field='hosts')
    aux = data[data['taxid1'] != data['taxid2']]
    common_prots = aux.drop_duplicates(["taxid1", "target"]).groupby(["target"],as_index=False).size()
    common_prots = common_prots[common_prots['size']>min_common]["target"].tolist()
    d = data[(data["source"].isin(common_prots)) | (data["target"].isin(common_prots))]
    build_graph(data=d, parasite=None, hosts=hosts, proteins=proteins, tissues=tissues, clusters=clusters, tags=tags, output_dir_path=output_dir_path)
    
def generate_gos_cytoscape_network(edges_file_path, functions, output_dir_path):
    with open(edges_file_path, 'r') as f:
        data = pd.read_csv(f, sep='\t')
    
    net = data[[0,1,2,3]]
    net = net[net['taxid1'] != net['taxid2']]
    
    go_net = []
    for ident in functions:
        gos = functions[ident]
        taxid = ident.split('.')[0]
        for go in gos:
            go_net.append((taxid, ident, 20, go))
    go_net = pd.DataFrame(go_net, columns=['taxid1', 'source', 'taxid2', 'target'])
    net = net.append(go_net)
    
    with open(os.path.join(output_dir_path, 'functional_network.tsv'), 'w') as outnet:
        net.to_csv(outnet, sep='\t', header=True, index=False, doublequote=False)
        
    
    
        
        
def build_graph(data, parasite, hosts, proteins, tissues, clusters, tags, output_dir_path):
    if parasite is not None:
        parasite_dir_path = os.path.join(output_dir_path, str(parasite))
        if not os.path.exists(parasite_dir_path):
            os.makedirs(parasite_dir_path)
        
        output_file_path = os.path.join(parasite_dir_path, 'dataset.json')
        d = data[data['taxid1'] == parasite]
        host_proteins = d['target'].tolist()
        d = data[(data['taxid1'] == parasite) | ((data['source'].isin(host_proteins)) & (data['target'].isin(host_proteins)))]
    else:
        output_file_path = os.path.join(output_dir_path, 'dataset.json')
        d = data.copy()
    net = None
    node_dict = {}
    edges = d.iloc[:, [1, 3, 6]].values.tolist() 
    for i, row in d.iterrows():
        source_tissues = []
        target_tissues = []
        if row["source"] in tissues:
            source_tissues = tissues[row["source"]]
        if row["target"] in tissues:
            target_tissues = tissues[row["target"]]
        tag1 = "Human protein" if row["taxid1"] in hosts else "Parasite protein"
        tag2 = "Human protein" if row["taxid2"] in hosts else "Parasite protein"
        source = {"key": row["source"],
                "label": proteins[row["source"]],
                "tag": tag1,
                "URL": "",
                "cluster": row["taxid1"],
                "tissues": source_tissues
                }
        target = {"key": row["target"],
                "label": proteins[row["target"]],
                "tag": tag2,
                "URL": "",
                "cluster": row["taxid2"],
                "tissues": target_tissues
                }
        node_dict.update({row['source']: source, row['target']: target})
    print("Building network")
    net = nx.Graph()
    net.add_weighted_edges_from(edges)
    nx.set_node_attributes(net, node_dict)
    print("Drawing layout")
    pos = nx.kamada_kawai_layout(net)
    print("Getting positions")
    coords = {}
    for n in pos:
        coords[n] = {"x": pos[n][0], "y": pos[n][1]}
    
    nx.set_node_attributes(net, coords)
    print("Calculating Centrality")
    centrality = nx.betweenness_centrality(net, weight='weight')
    nx.set_node_attributes(net, centrality, 'score')
    print("Saving graph")        
    graph = {"nodes": list(dict(net.nodes(data=True)).values()), "edges":edges, "clusters":clusters, "tags": tags}
    save_orthohpi_graph(graph=graph, output_file_path=output_file_path)
    print("saved")


def save_orthohpi_graph(graph, output_file_path):
    
    with open(output_file_path, 'w') as out:
        out.write(json.dumps(graph))
