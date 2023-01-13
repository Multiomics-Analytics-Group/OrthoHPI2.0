import os
import yaml
import json
import requests
import gzip
import itertools
import zipfile
import obonet
import networkx as nx
from Bio import SeqIO
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

def read_fasta(fasta_file_path):
    sequences = []
    fasta_sequences = SeqIO.parse(open(fasta_file_path),'fasta')
    for fasta in fasta_sequences:
        sequences.append(fasta.id)
    return sequences

def filter_sequences(sequences, valid_list):
    filter_out = []
    for parasite_id in valid_list:
        if parasite_id not in sequences:
            filter_out.append(parasite_id)
            
    return filter_out

def convert_df(df):
    return df.to_csv(sep='\t', header=True, index=False).encode('utf-8')

def export_graph(G, filename, format='graphml', output_dir='tmp'):
    file_path = os.path.join(output_dir, filename)
    if format == "graphml":
        nx.write_graphml_lxml(G, file_path)
    elif format == "cytoscape":
        cytoscape_data= nx.cytoscape_data(G)
        with open(file_path, 'w') as out:
            out.write(json.dumps(cytoscape_data))

def calculate_enrichment(pred_df, go_df):
    nodes = pred_df['source'].unique().tolist() + pred_df['target'].unique().tolist()
    total_nodes = len(nodes)
    selected_gos = go_df[go_df['#string_protein_id'].isin(nodes)].groupby('description').filter(lambda x: (len(x)> 10) & (len(x) < 500))['description'].unique().tolist()
    total_prots = len(go_df['#string_protein_id'].unique().tolist())
    enrichment = []
    for term in selected_gos:   
        members = go_df[(go_df['description'] == term)]['#string_protein_id']
        #E
        total_members = len(members)
        net_members = go_df[(go_df['description'] == term) & (go_df['#string_protein_id'].isin(nodes))]['#string_protein_id']
        #A
        total_net_members = len(net_members)
        
        odd_ratio, p_value = stats.fisher_exact([[total_net_members, total_nodes - total_net_members],
                                                [total_members - total_net_members, total_prots - total_members - total_nodes - total_net_members]])
        enrichment.append([term, total_net_members, total_nodes - total_net_members, total_members - total_net_members,  total_prots - total_members - total_nodes - total_net_members, p_value, odd_ratio, ','.join(net_members)])
    
    enrichment = pd.DataFrame(enrichment, columns=['go_term', 'A', 'B', 'C', 'D', 'p_value', 'odds_ratio', 'nodes'])
    if not enrichment.empty:
        enrichment['fdr_bh'] = multipletests(enrichment['p_value'].tolist(), alpha=0.01, method='fdr_bh')[1]
        enrichment = enrichment.sort_values(by='fdr_bh', ascending=True)
    
    return enrichment


def save_to_parquet(df, output_file):
    df.to_parquet(output_file, compression='gzip', index=False)


def read_parquet_file(input_file):
    df = pd.read_parquet(input_file)

    return df


def annotate_alias_id(predictions_df, taxids, config_file, sources, new_col, mapping_col):
    '''
    Adds an extra column to the provided dataframe with the String alias selected (e.g., UniProt id)

    :param DataFrame predictions_df: predictions dataframe to be annotated (requires mapping_col in columns)
    :param str config_file: path to config file (used to get the aliases for each species)
    :param list sources: what source ids need to be annotated

    :return DataFrame predictions_df: annotated dataframe with the String aliases of interest
    '''
    aliases = {}
    for taxid in taxids:
        aliases.update(parse_string_aliases(config_file=config_file, 
                    sources=sources, taxid=str(taxid), reverse=True))
    
    predictions_df[new_col] = predictions_df[mapping_col].map(aliases)
    #predictions_df['target_uniprot'] = predictions_df['target'].map(aliases)
    
    return predictions_df


def parse_string_aliases(config_file, sources, taxid='9606', reverse=False):
    '''
    Parses the alias file from String database and generates a dictionary
    that can be used to map to the right identifiers
    :param str config_file: path to the config file where the url to the String alias file should be defined
    :param list sources: list of sources that should be considered in the mapping (i.e. Ensembl_gene)
    :param str taxid: taxonomic identifier of the species for which to parse the aliases file
    :param bool reverse: whether to store alias --> string_id dictionary (False), or string_id --> alias (True)
    :return: dictionary with key --> alias, values --> string_id (reverse=False),
                or key --> string_id, values --> alias
    '''
    data_dict = {}
    urls = read_config(filepath=config_file, field='urls')
        
    if 'string_alias_url' in urls:
        filename = download_file(url=urls['string_alias_url'].replace('TAXID', taxid), data_dir='data')
    
    data = pd.read_csv(filename, sep='\t', header=0)
    if sources is not None:
        data = data[data['source'].isin(sources)]
    

    for i, row in data[['#string_protein_id', 'alias']].iterrows():
        if not reverse:
            data_dict[row['alias']] = row['#string_protein_id']
        else:
            data_dict[row['#string_protein_id']] = row['alias']
         
    return data_dict


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
    handle = gzip.open(filepath, "rb")

    return handle

def read_zipped_file(filepath):
    '''
    Opens a handler to access the content of zip file
    :param str filepath: path to the zip file
    :return: A bytes sequence that specifies the standard output
    '''
    file_name = filepath.split('/')[-1].split('.')[0]+'.tsv'
    archive = zipfile.ZipFile(filepath, 'r')
    handle = archive.open(file_name)

    return handle

def merge_dict_of_dicts(dict_of_dicts):
    dictionary = {}
    for d in dict_of_dicts:
        dictionary.update(dict_of_dicts[d])
        
    return dictionary
    
def merge_list_of_lists(list_of_lists):
    return list(itertools.chain.from_iterable(list_of_lists))


def convertOBOtoNet(ontologyFile):
    """
    Takes an .obo file and returns a NetworkX graph representation of the ontology, that holds multiple \
    edges between two nodes.
    :param str ontologyFile: path to ontology file.
    :return: NetworkX graph.
    """
    graph = obonet.read_obo(ontologyFile)

    return graph