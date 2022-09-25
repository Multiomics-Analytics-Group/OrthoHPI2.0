import os
import yaml
import requests
import gzip
import zipfile
from Bio import SeqIO
import pandas as pd


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


def parse_string_aliases(config_file, sources):
    '''
    Parses the alias file from String database and generates a dictionary
    that can be used to map to the right identifiers
    :param str config_file: path to the config file where the url to the String alias file should be defined
    :param list sources: list of sources that should be considered in the mapping (i.e. Ensembl_gene)
    :return: dictionary with key --> alias, values --> sources identifier
    '''
    data_dict = {}
    urls = read_config(filepath=config_file, field='urls')
    
    if 'string_alias_url' in urls:
        filename = download_file(url=urls['string_alias_url'], data_dir='data')
    
    data = pd.read_csv(filename, sep='\t', header=0)
    if sources is not None:
        data = data[data['source'].isin(sources)]
    

    for i, row in data[['#string_protein_id', 'alias']].iterrows():
        data_dict[row['alias']] = row['#string_protein_id']
         
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
    handle = gzip.open(filepath, "rt", encoding="utf8")

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
