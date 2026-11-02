'''Build mouse cell-type expression from the Tabula Muris Senis raw H5ADs.'''
import argparse
import os

import anndata
import numpy as np
import pandas as pd
from scipy import sparse

import utils

MOUSE_TAXID = '10090'
COUNTS_PER_MILLION = 1_000_000
CELLS_PER_CHUNK = 4_096
DEFAULT_AGE = '3m'
DOWNLOAD_DIR = os.path.join('data', 'downloads', 'tabula_muris_senis')
DEFAULT_INPUT = os.path.join(DOWNLOAD_DIR, 'tabula-muris-senis-droplet-official-raw-obj.h5ad')
DEFAULT_FACS_INPUT = os.path.join(DOWNLOAD_DIR, 'tabula-muris-senis-facs-official-raw-obj.h5ad')
# atlas tissue -> config['tissues'] labels. TISSUES records the large intestine under both
# `intestine` and `colon`
TISSUE_MAPPING = {
    'Bladder': ['urinary bladder'],
    'Brain': ['brain'],
    'Brain_Myeloid': ['brain'],
    'Brain_Non-Myeloid': ['brain'],
    'Heart': ['heart'],
    'Heart_and_Aorta': ['heart'],
    'Kidney': ['kidney'],
    'Large_Intestine': ['intestine', 'colon'],
    'Limb_Muscle': ['muscle'],
    'Liver': ['liver'],
    'Lung': ['lung'],
    'Marrow': ['bone marrow'],
    'Pancreas': ['pancreas'],
    'Skin': ['skin'],
    'Spleen': ['spleen'],
}


def source_tissues(obs):
    '''The atlas tissue of each cell, with Heart_and_Aorta split by its free annotation.'''
    source = obs['tissue'].astype(object).copy()
    if 'tissue_free_annotation' in obs:
        heart_or_aorta = source == 'Heart_and_Aorta'
        source.loc[heart_or_aorta] = (
            obs.loc[heart_or_aorta, 'tissue_free_annotation'].replace('', pd.NA)
            .fillna(source.loc[heart_or_aorta]).astype(object)
        )
    return source


def aggregate_expression(input_file, config_file, age):
    '''Pool the raw counts of each (tissue, cell type) and scale each pool to a million.'''
    atlas = anndata.read_h5ad(input_file, backed='r')
    obs = atlas.obs
    valid_tissues = {name.lower() for name in
                     utils.read_config(filepath=config_file, field='tissues').values()}

    # one row per (cell, label); a cell under two labels counts towards both
    cells = pd.DataFrame({
        'cell': np.arange(atlas.n_obs),
        'Tissue': source_tissues(obs).map(TISSUE_MAPPING).to_numpy(),
        'Cell type': obs['cell_ontology_class'].astype(object).to_numpy(),
    })
    cells = cells[(obs['age'].astype(str) == age).to_numpy()
                  & cells['Tissue'].notna().to_numpy() & cells['Cell type'].notna().to_numpy()]
    cells = cells.explode('Tissue')
    cells = cells[cells['Tissue'].isin(valid_tissues)]

    group_codes, groups = pd.factorize(pd.MultiIndex.from_frame(cells[['Tissue', 'Cell type']]))
    membership = sparse.csc_matrix(
        (np.ones(len(cells)), (group_codes, cells['cell'].to_numpy())),
        shape=(len(groups), atlas.n_obs))
    sums = np.zeros((len(groups), atlas.n_vars))
    for start in range(0, atlas.n_obs, CELLS_PER_CHUNK):
        stop = min(start + CELLS_PER_CHUNK, atlas.n_obs)
        block = membership[:, start:stop]
        if block.nnz:
            sums += (block @ atlas.X[start:stop].astype(np.float64)).toarray()

    aliases = utils.parse_string_aliases(
        config_file=config_file, sources=['Ensembl_gene', 'UniProt_GN_Name'], taxid=MOUSE_TAXID)
    genes = pd.Series(atlas.var_names.astype(str)).map(aliases)
    mapped = genes.notna().to_numpy()
    atlas.file.close()

    totals = sums.sum(axis=1, keepdims=True)
    cpm = np.divide(sums * COUNTS_PER_MILLION, totals, out=np.zeros_like(sums), where=totals > 0)
    data = pd.DataFrame(cpm[:, mapped], index=groups, columns=genes[mapped].to_numpy())
    data = data.rename_axis(['Tissue', 'Cell type']).reset_index()
    data = data.melt(id_vars=['Tissue', 'Cell type'], var_name='Gene', value_name='nCPM')
    data = data[data['nCPM'] > 0]
    # several atlas genes can map to one protein; keep the highest, as for HPA
    return (data.groupby(['Gene', 'Tissue', 'Cell type'], as_index=False)['nCPM'].max()
            [['Gene', 'Tissue', 'Cell type', 'nCPM']])


def combine_atlases(droplet, facs):
    '''Use FACS only for tissues absent from droplet data; never pool technologies.'''
    if facs is None:
        return droplet
    return pd.concat([droplet, facs[~facs['Tissue'].isin(droplet['Tissue'])]],
                     ignore_index=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', default=DEFAULT_INPUT, help='official raw droplet H5AD')
    parser.add_argument('--facs-input', default=DEFAULT_FACS_INPUT, help='official raw FACS H5AD')
    parser.add_argument('--droplet-only', action='store_true',
                        help='do not use FACS as a fallback for missing tissues')
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--output', default='data/mouse_atlas_cell_types.parquet')
    parser.add_argument('--age', default=DEFAULT_AGE, help='Tabula Muris Senis age (default: 3m)')
    args = parser.parse_args()

    droplet = aggregate_expression(input_file=args.input, config_file=args.config, age=args.age)
    facs = None
    if not args.droplet_only:
        facs = aggregate_expression(input_file=args.facs_input, config_file=args.config,
                                    age=args.age)
    data = combine_atlases(droplet, facs)
    utils.save_to_parquet(data, args.output)
    print(f'Wrote {len(data):,} mouse cell-type expression rows to {args.output}')
