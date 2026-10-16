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
import numpy as np
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

def prepare_go_background(go_df):
    '''
    What calculate_enrichment reads from a background, worked out once so a background
    tested against many protein sets is not regrouped for each of them.
    '''
    pairs = go_df[['term', '#string_protein_id']].drop_duplicates()

    return {'pairs': pairs,
            'annotated': set(pairs['#string_protein_id']),
            'names': dict(go_df[['term', 'description']].drop_duplicates().values),
            'sizes': pairs.groupby('term').size()}


def calculate_enrichment(proteins, go_df, min_term=10, max_term=500, max_share=0.25,
                         min_in_set=2, background=None):
    '''
    Fisher's exact test (one-sided, greater) of every Gene Ontology term against one set
    of proteins; `background` is prepare_go_background(go_df), when already worked out.
    '''
    if background is None:
        background = prepare_go_background(go_df)
    pairs, annotated = background['pairs'], background['annotated']
    sizes = background['sizes']
    # the universe is what is annotated
    tested = set(proteins) & annotated
    total_nodes = len(tested)
    total_prots = len(annotated)

    max_term = min(max_term, round(max_share * total_prots))

    hits = pairs[pairs['#string_protein_id'].isin(tested)]
    in_set = hits.groupby('term').size()
    terms = sizes[(sizes > min_term) & (sizes < max_term)].index.intersection(
        in_set[in_set >= min_in_set].index).sort_values()

    columns = ['go_id', 'go_term', 'A', 'B', 'C', 'D', 'p_value', 'odds_ratio', 'nodes']
    if len(terms) == 0:
        return pd.DataFrame(columns=columns)

    # 2x2 table: set proteins annotated to the term (A), rest of the set (B), outside
    # the set annotated (C), everything else (D)
    a = in_set[terms].to_numpy(dtype=np.int64)
    b = total_nodes - a
    c = sizes[terms].to_numpy(dtype=np.int64) - a
    d = total_prots - (a + c) - total_nodes + a
    # what scipy.stats.fisher_exact(alternative='greater') computes, over every term at once
    p_value = np.minimum(stats.hypergeom.cdf(b, a + b + c + d, a + b, b + d), 1.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        odds_ratio = np.where((c > 0) & (b > 0), a * d / (c * b), np.inf)
    empty = ((a + c) == 0) | ((b + d) == 0) | ((a + b) == 0) | ((c + d) == 0)
    p_value = np.where(empty, 1.0, p_value)
    odds_ratio = np.where(empty, np.nan, odds_ratio)

    members = (hits[hits['term'].isin(terms)].groupby('term')['#string_protein_id']
               .agg(lambda ids: ','.join(sorted(ids))))
    names = background['names']
    enrichment = pd.DataFrame({'go_id': terms,
                               'go_term': [names.get(t, t) for t in terms],
                               'A': a, 'B': b, 'C': c, 'D': d,
                               'p_value': p_value, 'odds_ratio': odds_ratio,
                               'nodes': members[terms].to_numpy()})
    # only the corrected p-values are used; the threshold is chosen in the app
    enrichment['fdr_bh'] = multipletests(enrichment['p_value'].tolist(), method='fdr_bh')[1]

    return enrichment.sort_values(by='fdr_bh', ascending=True)

def save_to_parquet(df, output_file):
    df.to_parquet(output_file, compression='gzip', index=False)


def read_parquet_file(input_file, filters=None):
    '''Reads a parquet file.'''
    df = pd.read_parquet(input_file, filters=filters)

    return df


def annotate_alias_id(predictions_df, taxids, config_file, sources, new_col, mapping_col):
    '''
    Adds an extra column to the provided dataframe with the String alias selected (e.g.,
    UniProt id)
    '''
    aliases = {}
    for taxid in taxids:
        aliases.update(parse_string_aliases(config_file=config_file, 
                    sources=sources, taxid=str(taxid), reverse=True))
    
    predictions_df[new_col] = predictions_df[mapping_col].map(aliases)
    
    return predictions_df


def parse_string_aliases(config_file, sources, taxid='9606', reverse=False):
    '''
    Parses the alias file from String database and generates a dictionary that can be
    used to map to the right identifiers
    '''
    data_dict = {}
    urls = read_config(filepath=config_file, field='urls')

    if 'string_alias_url' in urls:
        filename = download_file(url=urls['string_alias_url'].replace('TAXID', taxid), data_dir=os.path.join('data/downloads/species', str(taxid)))

    data = pd.read_csv(filename, sep='\t', header=0)
    if sources is not None:
        data = data[data['source'].isin(sources)]
        # later rows overwrite earlier ones, so descending preference leaves the most
        # preferred source written last; the sort is stable
        rank = {source: i for i, source in enumerate(sources)}
        data = data.sort_values('source', key=lambda s: s.map(rank), ascending=False, kind='stable')

    for i, row in data[['#string_protein_id', 'alias']].iterrows():
        if not reverse:
            data_dict[row['alias']] = row['#string_protein_id']
        else:
            data_dict[row['#string_protein_id']] = row['alias']
         
    return data_dict


def read_yaml(yaml_file):
    '''Reads YAML file and stores it in a dictionary'''
    content = None
    with open(yaml_file, 'r') as stream:
        try:
            content = yaml.safe_load(stream)
        except yaml.YAMLError as err:
            raise yaml.YAMLError("The yaml file {} could not be parsed. {}".format(yaml_file, err))
    return content

def read_config(filepath, field=None):
    '''
    Read the configuration file and return either the full content or an specific field.
    '''
    content = read_yaml(filepath)
    if content is not None:
        if field is not None:
            if field in content:
                return content[field]

    return content

def download_file(url, data_dir='data'):
    '''Download file from an url into an existing directory'''
    os.makedirs(data_dir, exist_ok=True)
    header = {'user-agent':'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/51.0.2704.103 Safari/537.36'}
    filename = url.split('/')[-1]
    filename = os.path.join(data_dir, filename)
    if not os.path.isfile(filename):
        r = requests.get(url, headers=header)
        with open(filename, 'wb') as out:
            out.write(r.content)
            
    return filename

def read_gzipped_file(filepath):
    '''
    Opens an underlying process to access a gzip file through the creation of a new pipe
    to the child.
    '''
    handle = gzip.open(filepath, "rb")

    return handle

def read_zipped_file(filepath):
    '''Opens a handler to access the content of zip file'''
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
    '''
    Takes an .obo file and returns a NetworkX graph representation of the ontology, that
    holds multiple     edges between two nodes.
    '''
    graph = obonet.read_obo(ontologyFile)

    return graph