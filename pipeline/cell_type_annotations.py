'''Load and combine cell-type annotations from host-specific atlas sources.'''

import pandas as pd
import utils

PIG_TAXID = '9823'
MOUSE_TAXID = '10090'
CELL_TYPE_COLUMNS = ['Gene', 'Tissue', 'Cell type', 'nTPM']

# HPA tissue names -> config['tissues'] labels. The gut tissues are written under their own
# label and `intestine`, since TISSUES annotates both; `pbmc` is HPA's only blood data (no
# erythrocytes or granulocytes)
HPA_TISSUE_LABELS = {
    'heart muscle': ['heart'],
    'skeletal muscle': ['muscle'],
    'bronchus': ['lung'],
    'pbmc': ['blood'],
    'colon': ['colon', 'intestine'],
    'rectum': ['rectum', 'intestine'],
    'small intestine': ['small intestine', 'intestine'],
}


def read_hpa(config_file):
    '''Reads the HPA file containing cell type protein expression profiles per tissue.'''
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
    '''Map gene identifiers and filetering only tissues relevant in OrthoHPI 2.0'''
    aliases = utils.parse_string_aliases(config_file, sources=['Ensembl_gene'])
    tissues = {t.lower() for t in utils.read_config(filepath=config_file, field='tissues').values()}
    # a tissue HPA reports keeps only the labels of config['tissues'] it stands for
    labels = {tissue: [label for label in HPA_TISSUE_LABELS.get(tissue, [tissue]) if label in tissues]
              for tissue in hpa_data['Tissue'].unique()}

    hpa_data = hpa_data.copy()
    hpa_data['Tissue'] = hpa_data['Tissue'].map(labels)
    hpa_data = hpa_data[hpa_data['Tissue'].str.len() > 0].explode('Tissue', ignore_index=True)
    hpa_data['Gene'] = hpa_data['Gene'].map(aliases)

    # several HPA tissues and Ensembl genes now fall under one label and protein; keep the
    # highest nTPM of each
    hpa_data = hpa_data.sort_values(by='nTPM', ascending=False).drop_duplicates(['Gene', 'Tissue', 'Cell type'], keep='first')

    return hpa_data


def filter_valid_proteins(hpa_data, valid_proteins):
    '''Keep only the HPA rows whose mapped Gene is in the pipeline's valid proteins.'''
    hpa_data = hpa_data[hpa_data['Gene'].isin(valid_proteins)]

    return hpa_data


def parse_hpa(config_file, valid_proteins):
    '''
    Load, map, and filter HPA single-cell data to (Gene, Tissue, Cell type, nTPM) rows
    for the given valid proteins.
    '''
    data = read_hpa(config_file=config_file)
    data = map_hpa_data(config_file=config_file, hpa_data=data)
    data = filter_valid_proteins(data, valid_proteins=valid_proteins)

    return data


def read_pig_atlas(data_dir, valid_proteins):
    '''Load preprocessed pig cell-atlas expression, when it has been generated.'''
    filename = f'{data_dir}/pig_atlas_cell_types.parquet'
    try:
        data = utils.read_parquet_file(input_file=filename)
    except FileNotFoundError:
        return pd.DataFrame(columns=CELL_TYPE_COLUMNS)

    missing = set(CELL_TYPE_COLUMNS).difference(data.columns)
    if missing:
        raise ValueError(f'{filename} is missing required columns: {sorted(missing)}')

    data = data[CELL_TYPE_COLUMNS]
    data = data[data['Gene'].astype(str).str.startswith(f'{PIG_TAXID}.')]
    return filter_valid_proteins(data, valid_proteins)


def read_mouse_atlas(data_dir, valid_proteins):
    '''Load preprocessed Tabula Muris Senis expression, when it has been generated.'''
    filename = f'{data_dir}/mouse_atlas_cell_types.parquet'
    try:
        data = utils.read_parquet_file(input_file=filename)
    except FileNotFoundError:
        return pd.DataFrame(columns=CELL_TYPE_COLUMNS)

    missing = set(CELL_TYPE_COLUMNS).difference(data.columns)
    if missing:
        raise ValueError(f'{filename} is missing required columns: {sorted(missing)}')

    data = data[CELL_TYPE_COLUMNS]
    data = data[data['Gene'].astype(str).str.startswith(f'{MOUSE_TAXID}.')]
    return filter_valid_proteins(data, valid_proteins)


def parse_cell_type_data(config_file, data_dir, valid_proteins):
    '''Return HPA and optional pig/mouse atlas cell-type annotations.'''
    human = parse_hpa(config_file=config_file, valid_proteins=valid_proteins)
    pig = read_pig_atlas(data_dir=data_dir, valid_proteins=valid_proteins)
    mouse = read_mouse_atlas(data_dir=data_dir, valid_proteins=valid_proteins)
    return pd.concat([human, pig, mouse], ignore_index=True, sort=False)
