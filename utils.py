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

def calculate_enrichment(proteins, go_df, min_term=10, max_term=500, min_in_set=2):
    """
    Fisher's exact test of every Gene Ontology term against one set of proteins.

    `proteins` and `go_df` have to be of the same species. The background is every
    annotated protein of `go_df`, so a set of one organism tested against the annotation of
    two is tested against a null that is partly another organism's: the two differ in how
    deeply they are annotated, and the difference is read as enrichment.

    A term is tested when the background gives it between `min_term` and `max_term`
    proteins -- the range that is neither nearly empty nor nearly everything -- and when at
    least `min_in_set` of the tested proteins carry it. The size range is read off the
    background rather than off the tested set: selecting a term because many of the tested
    proteins carry it and then testing whether they over-carry it is the same question
    asked twice, and on a set of twenty proteins it leaves almost nothing to test.

    :param proteins: the protein ids to test, as STRING ids
    :param go_df: GO annotations of their species, with columns #string_protein_id and
                  description
    :param int min_term: smallest background term tested, exclusive
    :param int max_term: largest background term tested, exclusive
    :param int min_in_set: proteins of the set a term needs before it is worth testing
    :return: dataframe of go_term, the four cells of the table, p_value, odds_ratio, the
             proteins behind it, and the Benjamini-Hochberg corrected fdr_bh
    """
    annotated = set(go_df['#string_protein_id'])
    # the universe is what is annotated: a protein no term can be counted against belongs
    # in neither margin of the table
    tested = set(proteins) & annotated
    total_nodes = len(tested)
    total_prots = len(annotated)

    sizes = go_df.groupby('description')['#string_protein_id'].nunique()
    in_set = (go_df[go_df['#string_protein_id'].isin(tested)]
              .groupby('description')['#string_protein_id'].nunique())
    terms = sizes[(sizes > min_term) & (sizes < max_term)].index.intersection(
        in_set[in_set >= min_in_set].index)

    enrichment = []
    for term, ids in go_df[go_df['description'].isin(terms)].groupby(
            'description')['#string_protein_id']:
        members = set(ids)
        net_members = members & tested

        # the 2x2 table the test is run on: the proteins of the set annotated to the term
        # (A), the rest of the set (B), the proteins outside the set annotated to it (C)
        # and everything else (D). The set's proteins annotated to the term are taken out
        # of both margins of D, so they are added back once
        a = len(net_members)
        b = total_nodes - a
        c = len(members) - a
        d = total_prots - len(members) - total_nodes + a
        odd_ratio, p_value = stats.fisher_exact([[a, b], [c, d]])
        enrichment.append([term, a, b, c, d, p_value, odd_ratio, ','.join(sorted(net_members))])

    enrichment = pd.DataFrame(enrichment, columns=['go_term', 'A', 'B', 'C', 'D', 'p_value',
                                                   'odds_ratio', 'nodes'])
    if not enrichment.empty:
        enrichment['fdr_bh'] = multipletests(enrichment['p_value'].tolist(), alpha=0.01,
                                             method='fdr_bh')[1]
        enrichment = enrichment.sort_values(by='fdr_bh', ascending=True)

    return enrichment

def save_to_parquet(df, output_file):
    df.to_parquet(output_file, compression='gzip', index=False)


def read_parquet_file(input_file, filters=None):
    '''
    Reads a parquet file. Optional filters (i.e. [('taxid', 'in', [9606])]) are
    pushed down to the reader so only the matching row groups are loaded.
    '''
    df = pd.read_parquet(input_file, filters=filters)

    return df


def annotate_alias_id(predictions_df, taxids, config_file, sources, new_col, mapping_col):
    '''
    Adds an extra column to the provided dataframe with the String alias selected (e.g., UniProt id)

    :param DataFrame predictions_df: predictions dataframe to be annotated (requires mapping_col in columns)
    :param str config_file: path to config file (used to get the aliases for each species)
    :param list sources: what source ids need to be annotated, in order of preference
                (see parse_string_aliases)

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
    :param list sources: sources to consider (i.e. Ensembl_gene), in order of preference:
                the alias of the first source that has one for an identifier wins. STRING
                alias files differ per species -- only human carries
                Ensembl_HGNC_uniprot_ids -- so listing fallbacks maps the other species
                without changing the ones that do have the preferred source.
    :param str taxid: taxonomic identifier of the species for which to parse the aliases file
    :param bool reverse: whether to store alias --> string_id dictionary (False), or string_id --> alias (True)
    :return: dictionary with key --> alias, values --> string_id (reverse=False),
                or key --> string_id, values --> alias
    '''
    data_dict = {}
    urls = read_config(filepath=config_file, field='urls')

    if 'string_alias_url' in urls:
        filename = download_file(url=urls['string_alias_url'].replace('TAXID', taxid), data_dir=os.path.join('data/downloads/species', str(taxid)))

    data = pd.read_csv(filename, sep='\t', header=0)
    if sources is not None:
        data = data[data['source'].isin(sources)]
        # Rows are written into data_dict in order and later ones overwrite earlier
        # ones, so sorting by descending preference leaves the most preferred source
        # written last -- and therefore the one that wins. The sort is stable, so
        # within a single source the file order still decides, as it did before.
        rank = {source: i for i, source in enumerate(sources)}
        data = data.sort_values('source', key=lambda s: s.map(rank), ascending=False, kind='stable')

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