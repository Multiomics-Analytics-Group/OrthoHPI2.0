import pandas as pd
import utils


def read_hpa(config_file):
    '''
    Reads the HPA file containing cell type protein expression profiles
    per tissue.
    :param str config_file: path to the configuration file
    :return: pandas dataframe with the protein expression profiles for each tissue and cell type
    '''
    urls = utils.read_config(filepath=config_file, field='urls')
    if 'hpa_single_cell_tissue_url' not in urls:
        raise KeyError("hpa_single_cell_tissue_url missing from config urls")

    filename = utils.download_file(url=urls['hpa_single_cell_tissue_url'], data_dir='data/downloads')
    hpa_file = utils.read_zipped_file(filepath=filename)
    data = pd.read_csv(hpa_file, sep='\t', header=0)
    data = data.sort_values(by='nTPM', ascending=False).drop_duplicates(['Gene', 'Tissue', 'Cell type'], keep='first')
    data = data[data['nTPM'] > 0.0]

    return data


def map_hpa_data(config_file, hpa_data):
    '''
    Map gene identifiers and filetering only tissues relevant in OrthoHPI 2.0
    :param str config_file: path to the config file
    :param dataframe hpa_data: pandas dataframe with the single cell type data from HPA
    :return: mapped dataframe
    '''
    aliases = utils.parse_string_aliases(config_file, sources=['Ensembl_gene'])
    tissues_mapping = {'heart muscle':'heart', 'small intestine':'intestine', 'rectum':'intestine', 'bronchus':'lung', 'colon':'intestine'}
    hpa_data = hpa_data.copy()
    hpa_data['Tissue'] = hpa_data['Tissue'].replace(tissues_mapping)
    hpa_data['Gene'] = hpa_data['Gene'].map(aliases)

    tissues = utils.read_config(filepath=config_file, field='tissues')
    hpa_data = hpa_data[hpa_data['Tissue'].isin([t.lower() for t in tissues.values()])]

    return hpa_data


def filter_valid_proteins(hpa_data, valid_proteins):
    '''
    Keep only the HPA rows whose mapped Gene is in the pipeline's valid proteins.
    :param dataframe hpa_data: mapped HPA dataframe (Gene column holds STRING protein ids)
    :param iterable valid_proteins: protein ids kept after the pipeline filters
    :return: filtered dataframe
    '''
    hpa_data = hpa_data[hpa_data['Gene'].isin(valid_proteins)]

    return hpa_data


def parse_hpa(config_file, valid_proteins):
    '''
    Load, map, and filter HPA single-cell data to (Gene, Tissue, Cell type, nTPM)
    rows for the given valid proteins.
    :param str config_file: path to the configuration file
    :param iterable valid_proteins: protein ids kept after the pipeline filters
    :return: filtered HPA dataframe
    '''
    data = read_hpa(config_file=config_file)
    data = map_hpa_data(config_file=config_file, hpa_data=data)
    data = filter_valid_proteins(data, valid_proteins=valid_proteins)

    return data
