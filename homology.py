from email.errors import FirstHeaderLineIsContinuationDefect
import os
import yaml
import requests
import gzip
import pandas as pd
import networkx as nx
import json
from joblib import Parallel, delayed

def read_yaml(yaml_file):
    """
    Reads YAML file and stores it in a dictionary
    :param str yaml_file: path to yaml file
    :return: a dictionary with the content of the yaml file
    """
    content = None
    with open(yaml_file, 'r') as stream:
        try:
            content = yaml.safe_load(stream)
        except yaml.YAMLError as err:
            raise yaml.YAMLError("The yaml file {} could not be parsed. {}".format(yaml_file, err))
    return content

def read_config(filepath, field=None):
    """
    Read the configuration file and return either the full content or an specific field.
    
    :param str filepath: path to configuration file
    :param str field: field to be obtained from the configuration
    
    :return: dictionary with the content of the configuration or the field specified
    """
    content = read_yaml(filepath)
    if content is not None:
        if field is not None:
            if field in content:
                return content[field]

    return content

def download_file(url, data_dir='data'):
    """
    Download file from an url into an existing directory
    :param str url: URL address where to download the data from
    :param str data_dir: path to directory where to download the data
    :return: filepath to the downloaded data
    """
    header = {'user-agent':'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/51.0.2704.103 Safari/537.36'}
    filename = url.split('/')[-1]
    filename = os.path.join(data_dir, filename)
    if not os.path.isfile(filename):
        r = requests.get(url, headers=header)
        with open(filename, 'wb') as out:
            out.write(r.content)
            
    return filename

def read_gzipped_file(filepath):
    """
    Opens an underlying process to access a gzip file through the creation of a new pipe to the child.
    :param str filepath: path to gzip file.
    :return: A bytes sequence that specifies the standard output.
    """
    handle = gzip.open(filepath, "rt")

    return handle

def merge_dict_of_dicts(dict_of_dicts):
    dictionary = {}
    for d in dict_of_dicts:
        dictionary.update(dict_of_dicts[d])
        
    return dictionary
    

def get_eggnog_orthologs(url, proteins, taxids, group):
    """
    NOT USED
    Obtains the list of orthologs and paralogs for a list of proteins and taxids for a specific eggNOG group
    :param str url: EggNOG API address where to retrieve fine grained orthologs and paralogs response ----> [[[taxid, species, [paralogs], [[taxid1, species1, [orthologs1]], [[taxid2, species2, [orthologs2]]]]]
    :param list proteins: list of proteins for which to obtain orthologs and paralogs for the provided taxids
    :param list taxids: list of taxonomic identifiers for which to obtain orthologs and paralogs
    :param str group: EggNOG group for which to obtain orthologs and paralogs
    :return: two lists with orthologs and paralogs identified
    """
    orthologs = []
    paralogs = []
    url = url.replace('PROTEINS', ",".join(proteins)).replace('TAXIDS',",".join(taxids)).replace('GROUP', group)
    response = requests.get(url)
    data = response.json()
    if 'orthologs' in data:
        if len(data['orthologs']) > 0:
            paralogs = data['orthologs'][0][0][2]
            orthologs = data['orthologs'][0][1]

    return orthologs, paralogs

def get_proteins(config_file):
    """
    Retrieve all proteins for all species
    :param str config_file: path to config file
    
    :return: dictionary with all proteins for all species. Key -> tax id, value -> dictionary: key -> protein id, value -> protein name"""
    proteins = {}
    hosts = read_config(filepath=config_file, field='hosts')
    parasites = read_config(filepath=config_file, field='parasites')
    urls = read_config(filepath=config_file, field='urls')
    if "string_protein_url" in urls:
        string_file = urls['string_protein_url']
        if hosts is not None and parasites is not None:
            taxids = list(hosts.keys()) + list(parasites.keys())
            for taxid in taxids:
                proteins[taxid] = parse_proteins(string_file, taxid)
    
    return proteins

def parse_proteins(string_file, taxid):
    """
    Retrieve proteins for a given specie
    :param str string_file: url to string PPI file
    :param str config_file: path to config file
    
    :return: dictionary with all proteins. Key -> Ensembl protein id, value -> protein name
    """
    proteins = {}
    if string_file is not None:
        filename = download_file(url=string_file.replace('TAXID', str(taxid)), data_dir='data')
        sp = read_gzipped_file(filename)
        first = True
        for line in sp:
            if first:
                first = False
                continue
            
            data = line.rstrip().split('\t')
            identifier = data[0]
            name = data[1]
            proteins[identifier] = name
            
    return proteins

def get_host_ppi(config_file, valid_proteins, ouput_filepath):
    """
    Retrieve Host intra-species proteins
    :param str config_file: path to config file
    :param dict valid_proteins: dictionary with valid proteins
    :param str ouput_filepath: path to edges file
    protein1 protein2 neighborhood fusion cooccurence coexpression experimental database textmining combined_score
    """
    seen = set()
    hosts = read_config(filepath=config_file, field='hosts')
    urls = read_config(filepath=config_file, field='urls')
    with open(ouput_filepath, 'a') as out:
        if "string_ppi_url" in urls:
            string_file = urls['string_ppi_url']
            if hosts is not None:
                for host in hosts:
                    if string_file is not None:
                        filename = download_file(url=string_file.replace('HOST', str(host)), data_dir='data')
                        sp = read_gzipped_file(filename)
                        first = True
                        for line in sp:
                            if first:
                                first = False
                                continue
                            data = line.rstrip().split(' ')
                            protein1 = data[0]
                            protein2 = data[1]
                            exp_score = float(data[6])/1000
                            db_score = float(data[7])/1000
                            avg_score = (exp_score + db_score)/2
                            if protein1 in valid_proteins and protein2 in valid_proteins:
                                if exp_score >= 0.7 or db_score >= 0.7:
                                    if (protein1, protein2) not in seen:
                                        out.write("\t".join([str(host), protein1, str(host), protein2, str(exp_score), str(db_score), str(avg_score), "", "", "intra-species"])+"\n")
                                        seen.add((protein1, protein2))
                                        seen.add((protein2, protein1))

        
def get_eggnog_groups(filepath, proteins):
    """
    Obtains all the EggNOG groups which contains a list of given proteins
    :param str filepath: path to the EggNOG groups file
    :param list proteins: list of Ensembl protein identifiers
    :return: dictionary with all the valid EggNOG groups. Key -> group, value -> list proteins in the group
    """
    valid_groups = {}
    groups = read_gzipped_file(filepath)
    first = True
    for line in groups:
        if first:
            first = False
            continue
        data = line.rstrip().split('\t')
        group = data[1]
        gproteins = data[4].split(',')
        int_proteins = list(set(proteins).intersection(gproteins))
        if  len(int_proteins) > 0:
            valid_groups[group] = int_proteins
    
    return valid_groups
            
def get_links(filepath, valid_groups, ouput_filepath, config_file):
    """
    Obtain the transferred interactions at the EggNOG group level from STRING
    Writes into a file 'predictions.tsv' with the list of predicted links based on homology.
    Structure of the file: taxid1, protein1, taxid2, protein2, str(experimental_evidence), str(databases_evidence), group1, group2
    :param str filepath: path to STRING file with the groups links
    :param dict valid_groups: dictionary with all the valid groups
    :param str output_filepath: path to output file
    """
    links = []
    seen = set()
    cog_links = read_gzipped_file(filepath)
    hosts = read_config(filepath=config_file, field='hosts')
    first = True
    with open(ouput_filepath, 'w') as out:
        out.write("\t".join(["taxid1", "source", "taxid2", "target", "experimental_evidence_score", "databases_evidence_score", "average_score", "group1", "group2", "edge_type"])+"\n")
        for line in cog_links:
            if first:
                first = False
                continue
            data = line.rstrip().split(' ')
            group1 = data[0]
            group2 = data[1]
            experimental_evidence = int(data[6])/1000
            databases_evidence = int(data[7])/1000
            
            if group1 in valid_groups and group2 in valid_groups:
                if experimental_evidence >= 0.7 or databases_evidence >= 0.7:
                    average_score = (experimental_evidence + databases_evidence) / 2
                    for protein1 in valid_groups[group1]:
                        for protein2 in valid_groups[group2]:
                            taxid1 = protein1.split('.')[0]
                            taxid2 = protein2.split('.')[0]
                            if (int(taxid1) in hosts or int(taxid2) in hosts) and (taxid1 != taxid2):
                                if (protein1, protein2) not in seen:
                                    if int(taxid1) in hosts:
                                        aux_id = taxid1
                                        aux_prot = protein1
                                        aux_group = group1
                                        taxid1 = taxid2
                                        protein1 = protein2
                                        group1 = group2
                                        taxid2 = aux_id
                                        protein2 = aux_prot
                                        group2 = aux_group
                                    out.write("\t".join([taxid1, protein1, taxid2, protein2, str(experimental_evidence), str(databases_evidence), str(average_score), group1, group2, "inter-species"])+"\n")
                                    seen.add((protein1, protein2))
                                    seen.add((protein2, protein1))
                                    
def setup(config_file, output_file_path):
    """
    Downloads all necessary files according to the urls specified in configuration
    except the protein files and host PPIs
    
    :param str config_file: path to the configuration file
    """
    urls = read_config(filepath=config_file, field='urls')
    for url_name in urls:
        url = urls[url_name]
        if url_name != "string_protein_url" and url_name != "string_ppi_url":
            filename = download_file(url=url, data_dir=output_file_path)


def apply_tissue_filter(config_file, valid_proteins, cutoff):
    hosts = read_config(filepath=config_file, field='hosts')
    for taxid in hosts:
        proteins = valid_proteins[taxid]
        #print("T before", len(proteins))
        if 'tissues_url' in hosts[taxid]:
            url = hosts[taxid]['tissues_url']
            filename = download_file(url=url, data_dir='data')
            tissues, proteins = get_tissues(filename, proteins, cutoff)
    
        valid_proteins[taxid] = proteins
        #print("T after", len(valid_proteins[taxid]))

    return tissues

def get_tissues(tissues_file, valid_proteins, cutoff):
    """
    Get protein tissue expression for relevant tissues in the lifecycle of the
    studied parasites
    
    :param str tissues_file: path to file with tissue expression (tissues.jensenlab.org)
    :param list valid_proteins: all proteins studied
    :param float cutoff: minimum confidence score accepted (tissues.jensenlab.org)
    
    :return tissues: dictionary with protein tissue expression
    """
    tissues = {}
    filters = {}
    valid_tissues = set()
    parasites = read_config(filepath=config_file, field='parasites')
    for parasite in parasites:
        t = parasites[parasite]['tissues']
        valid_tissues.update(t)
    
    first = True
    with open(tissues_file, 'r') as f:
        for line in f:
            if first:
                first = False
                continue
            
            data = line.rstrip().split('\t')
            protein  = "9606."+data[0]
            tissue = data[2]
            score = float(data[6])
            
            if protein in valid_proteins and score > cutoff and tissue in valid_tissues:
                if protein not in tissues:
                    tissues[protein] = []
                
                tissues[protein].append(tissue)
                filters[protein] = valid_proteins[protein]
    
    return tissues, filters


def apply_compartment_filter(config_file, valid_proteins, cutoff):
    hosts = read_config(filepath=config_file, field='hosts')
    for taxid in hosts:
        proteins = valid_proteins[taxid]
        #print("C before", len(proteins))
        if 'compartments_url' in hosts[taxid]:
            url = hosts[taxid]['compartments_url']
            filename = download_file(url=url, data_dir='data')
            compartments, proteins = get_compartments(filename, proteins, cutoff)
    
        valid_proteins[taxid] = proteins
        #print("C after", len(valid_proteins[taxid]))
    
    return compartments


def get_compartments(compartments_file, valid_proteins, cutoff):
    """
    Get protein cellular compartment expression relevant in the lifecycle of the
    studied parasites
    
    :param str compartments_file: path to file with cellular compartment expression (compartments.jensenlab.org)
    :param dict valid_proteins: dictionary with annotations in valid proteins
    :param float cutoff: minimum confidence score accepted (tissues.jensenlab.org)
    
    :return filtered_dict: dictionary with only proteins in relevant compartments
    """
    compartments = {}
    filters = {}
    valid_compartments = set()
    parasites = read_config(filepath=config_file, field='parasites')
    for parasite in parasites:
        if "compartments" in parasites[parasite]:
            t = parasites[parasite]['compartments']
        else:
            t = 'GO:0005886'
        valid_compartments.add(t)
    
    first = True
    with open(compartments_file, 'r') as f:
        for line in f:
            if first:
                first = False
                continue
            
            data = line.rstrip().split('\t')
            protein  = "9606."+data[0]
            compartment = data[2]
            score = float(data[4])
            if protein in valid_proteins and score > cutoff and compartment in valid_compartments:
                if protein not in compartments:
                    compartments[protein] = []
                
                compartments[protein].append(compartment)
                filters[protein] = valid_proteins[protein]
  
    return compartments, filters


def apply_fileters(proteins, annotations):
    """
    Filter out proteins"""

def get_clusters_from_config(config_file):
    """"clusters": [{ "key": string, "color": string, "clusterLabel": string }], #key:tax id, color:hex, clusterLabel: name"""
    clusters = []
    hosts = read_config(filepath=config_file, field='hosts')
    parasites = read_config(filepath=config_file, field='parasites')
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
        
    parasites = read_config(filepath=config_file, field='parasites')
    hosts = read_config(filepath=config_file, field='hosts')
    Parallel(n_jobs=n)(delayed(build_graph)(data, parasite, hosts, proteins, tissues, clusters, tags, output_dir_path) for parasite in parasites)
    
    
def generate_common_cytoscape_network(edges_file_path, proteins, tissues, config_file, output_dir_path, min_common=2):
    clusters = get_clusters_from_config(config_file)
    tags = [{"key": "Parasite protein", "image": "concept.svg"}, {"key":"Human protein", "image": "person.svg" }]
    with open(edges_file_path, 'r') as f:
        data = pd.read_csv(f, sep='\t')
        
    hosts = read_config(filepath=config_file, field='hosts')
    aux = data[data['taxid1'] != data['taxid2']]
    common_prots = aux.drop_duplicates(["taxid1", "target"]).groupby(["target"],as_index=False).size()
    common_prots = common_prots[common_prots['size']>min_common]["target"].tolist()
    d = data[(data["source"].isin(common_prots)) | (data["target"].isin(common_prots))]
    build_graph(data=d, parasite=None, hosts=hosts, proteins=proteins, tissues=tissues, clusters=clusters, tags=tags, output_dir_path=output_dir_path)
    
        
        
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


if __name__ == "__main__":
    data_dir = 'data'
    config_file = 'config.yml'
    
    setup(config_file=config_file, output_file_path=data_dir)
    proteins = get_proteins(config_file)
    tissues = apply_tissue_filter(config_file, proteins, cutoff=3.0)
    compartments = apply_compartment_filter(config_file, proteins, cutoff=3.5)
    proteins = merge_dict_of_dicts(dict_of_dicts=proteins)
    valid_groups = get_eggnog_groups(filepath=os.path.join(data_dir, '2759_members.tsv.gz'), proteins=proteins.keys())
    get_links(filepath=os.path.join(data_dir, 'COG.links.detailed.v11.5.txt.gz'), valid_groups=valid_groups,
              ouput_filepath=os.path.join(data_dir, 'predictions.tsv'), config_file=config_file)
    get_host_ppi(config_file=config_file, valid_proteins=proteins, ouput_filepath=os.path.join(data_dir, 'predictions.tsv'))
    generate_cytoscape_network(edges_file_path=os.path.join(data_dir, 'predictions.tsv'), proteins=proteins,
                                        tissues=tissues, config_file=config_file, output_dir_path='web/public', n=6)
    generate_common_cytoscape_network(edges_file_path=os.path.join(data_dir, 'predictions.tsv'), proteins=proteins,
                                        tissues=tissues, config_file=config_file, output_dir_path='web/public', min_common=15)
    
    
    
