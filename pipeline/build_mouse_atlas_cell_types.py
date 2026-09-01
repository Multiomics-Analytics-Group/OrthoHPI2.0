"""
Build mouse cell-type expression from the Tabula Muris Senis droplet H5AD.

By default this downloads the official raw droplet dataset from Figshare, filters
to 3-month mice, normalizes each cell to 10,000 counts, applies log1p, averages
expression by tissue and Cell Ontology label, and maps mouse genes to STRING IDs:

    .venv/bin/python -m pipeline.build_mouse_atlas_cell_types

Supply --input to process a previously downloaded H5AD. The default droplet
dataset is used deliberately: its UMI counts must not be pooled directly with
the Smart-seq2/FACS dataset's read counts.
"""
import argparse
import os

import anndata
import h5py
import numpy as np
import pandas as pd
import requests
from scipy import sparse

import utils

MOUSE_TAXID = '10090'
COUNTS_PER_CELL = 10_000
CELLS_PER_CHUNK = 4_096
DEFAULT_AGE = '3m'
DEFAULT_INPUT = os.path.join('data', 'downloads', 'tabula_muris_senis',
                             'tabula-muris-senis-droplet-official-raw-obj.h5ad')
DEFAULT_FACS_INPUT = os.path.join('data', 'downloads', 'tabula_muris_senis',
                                  'tabula-muris-senis-facs-official-raw-obj.h5ad')
HDF5_SIGNATURE = b'\x89HDF\r\n\x1a\n'
TISSUE_MAPPING = {
    'Bladder': 'urinary bladder',
    'Brain': 'brain',
    'Brain_Myeloid': 'brain',
    'Brain_Non-Myeloid': 'brain',
    'Heart': 'heart',
    'Heart_and_Aorta': 'heart',
    'Kidney': 'kidney',
    'Large_Intestine': 'intestine',
    'Limb_Muscle': 'muscle',
    'Liver': 'liver',
    'Lung': 'lung',
    'Marrow': 'bone marrow',
    'Pancreas': 'pancreas',
    'Skin': 'skin',
    'Spleen': 'spleen',
}


def normalized_tissues(obs):
    """Map source tissue labels to the OrthoHPI vocabulary without inventing matches."""
    source = obs['tissue'].astype(object).copy()
    if 'tissue_free_annotation' in obs:
        # Only this pooled source category needs its free annotation to distinguish
        # heart from aorta. Other free annotations are finer tissue names that are
        # not part of the OrthoHPI tissue vocabulary.
        heart_or_aorta = source == 'Heart_and_Aorta'
        source.loc[heart_or_aorta] = (
            obs.loc[heart_or_aorta, 'tissue_free_annotation'].replace('', pd.NA)
            .fillna(source.loc[heart_or_aorta]).astype(object)
        )
    return source.map(TISSUE_MAPPING)


def is_hdf5(filename):
    """Return whether filename starts with the HDF5 file signature."""
    try:
        with open(filename, 'rb') as handle:
            return handle.read(len(HDF5_SIGNATURE)) == HDF5_SIGNATURE
    except FileNotFoundError:
        return False


def download_atlas(url, output_file):
    """Stream the large atlas download into a complete, validated H5AD file."""
    if is_hdf5(output_file):
        return output_file

    if os.path.exists(output_file):
        print(f'Removing invalid atlas download: {output_file}')
        os.remove(output_file)
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    temporary_file = f'{output_file}.part'
    if os.path.exists(temporary_file):
        os.remove(temporary_file)

    print(f'Downloading Tabula Muris Senis atlas to {output_file}...')
    complete = False
    try:
        with requests.get(url, stream=True, timeout=60) as response:
            response.raise_for_status()
            with open(temporary_file, 'wb') as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
        if not is_hdf5(temporary_file):
            raise ValueError('Downloaded file is not a valid HDF5 file')
        os.replace(temporary_file, output_file)
        complete = True
    finally:
        if not complete and os.path.exists(temporary_file):
            os.remove(temporary_file)
    return output_file


def mean_normalized_expression(matrix):
    """Return log1p(CP10K)-normalized mean expression for every gene."""
    matrix = matrix.tocsr().astype(np.float64)
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    nonzero = totals > 0
    matrix = matrix[nonzero]
    if matrix.shape[0] == 0:
        return np.zeros(matrix.shape[1])

    matrix = sparse.diags(COUNTS_PER_CELL / totals[nonzero]) @ matrix
    matrix.data = np.log1p(matrix.data)
    return np.asarray(matrix.mean(axis=0)).ravel()


def read_csr_chunk(matrix, start, stop, n_genes):
    """Read consecutive rows from an H5AD CSR matrix into an in-memory CSR matrix."""
    offsets = matrix['indptr'][start:stop + 1]
    data_start, data_stop = offsets[0], offsets[-1]
    return sparse.csr_matrix(
        (matrix['data'][data_start:data_stop], matrix['indices'][data_start:data_stop],
         offsets - data_start),
        shape=(stop - start, n_genes),
    )


def normalized_expression_sum(matrix):
    """Return summed log1p(CP10K) expression and the number of nonempty cells."""
    matrix = matrix.tocsr().astype(np.float64)
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    nonzero = totals > 0
    matrix = matrix[nonzero]
    if matrix.shape[0] == 0:
        return np.zeros(matrix.shape[1]), 0

    matrix = sparse.diags(COUNTS_PER_CELL / totals[nonzero]) @ matrix
    matrix.data = np.log1p(matrix.data)
    return np.asarray(matrix.sum(axis=0)).ravel(), matrix.shape[0]


def aggregate_expression(input_file, config_file, age):
    """Aggregate the selected Tabula Muris Senis H5AD into the shared schema."""
    atlas = anndata.read_h5ad(input_file, backed='r')
    required = {'age', 'tissue', 'cell_ontology_class'}
    missing = required.difference(atlas.obs.columns)
    if missing:
        raise ValueError(f'{input_file} is missing required obs columns: {sorted(missing)}')

    obs_columns = list(required | {'tissue_free_annotation'})
    obs_columns = [column for column in obs_columns if column in atlas.obs.columns]
    obs = atlas.obs[obs_columns].copy()
    obs = obs[obs['age'].astype(str) == age]
    obs['Tissue'] = normalized_tissues(obs)
    valid_tissues = {name.lower() for name in
                     utils.read_config(filepath=config_file, field='tissues').values()}
    obs = obs[obs['Tissue'].isin(valid_tissues) & obs['cell_ontology_class'].notna()]

    aliases = utils.parse_string_aliases(
        config_file=config_file,
        sources=['Ensembl_gene', 'UniProt_GN_Name'],
        taxid=MOUSE_TAXID,
    )
    genes = pd.Series(atlas.var_names.astype(str), index=np.arange(atlas.n_vars)).map(aliases)
    mapped = genes.notna().to_numpy()

    group_codes, groups = pd.factorize(pd.MultiIndex.from_frame(
        obs[['Tissue', 'cell_ontology_class']]))
    cell_groups = np.full(atlas.n_obs, -1, dtype=np.int32)
    cell_groups[atlas.obs.index.get_indexer(obs.index)] = group_codes
    sums = [np.zeros(atlas.n_vars) for _ in groups]
    counts = np.zeros(len(groups), dtype=np.int64)
    n_obs = atlas.n_obs
    n_vars = atlas.n_vars
    atlas.file.close()

    with h5py.File(input_file, 'r') as h5ad:
        matrix = h5ad['X']
        if matrix.attrs.get('h5sparse_format') != 'csr':
            raise ValueError(f'{input_file} must store X as a CSR sparse matrix')
        for start in range(0, n_obs, CELLS_PER_CHUNK):
            stop = min(start + CELLS_PER_CHUNK, n_obs)
            chunk_groups = cell_groups[start:stop]
            if not (chunk_groups >= 0).any():
                continue
            chunk = read_csr_chunk(matrix, start, stop, n_vars)
            for group_code in np.unique(chunk_groups[chunk_groups >= 0]):
                expression, cell_count = normalized_expression_sum(
                    chunk[chunk_groups == group_code])
                sums[group_code] += expression
                counts[group_code] += cell_count

    rows = []
    for group_code, (tissue, cell_type) in enumerate(groups):
        if counts[group_code] == 0:
            continue
        means = (sums[group_code] / counts[group_code])[mapped]
        frame = pd.DataFrame({
            'Gene': genes[mapped].to_numpy(),
            'Tissue': tissue,
            'Cell type': str(cell_type),
            'nTPM': means,
        })
        rows.append(frame[frame['nTPM'] > 0])

    if not rows:
        return pd.DataFrame(columns=['Gene', 'Tissue', 'Cell type', 'nTPM'])
    data = pd.concat(rows, ignore_index=True)
    return (data.groupby(['Gene', 'Tissue', 'Cell type'], as_index=False, observed=True)['nTPM']
            .mean())


def build(input_file, config_file, output_file, age):
    """Write Tabula Muris Senis cell-type expression as a shared Parquet artifact."""
    data = aggregate_expression(input_file=input_file, config_file=config_file, age=age)
    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    utils.save_to_parquet(data, output_file)
    print(f'Wrote {len(data):,} mouse cell-type expression rows to {output_file}')
    return data


def combine_atlases(droplet, facs):
    """Use FACS only for tissues absent from droplet data; never pool technologies."""
    if facs is None:
        return droplet
    return pd.concat([droplet, facs[~facs['Tissue'].isin(droplet['Tissue'])]],
                     ignore_index=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', default=None, help='official raw droplet H5AD')
    parser.add_argument('--facs-input', default=None, help='official raw FACS H5AD')
    parser.add_argument('--droplet-only', action='store_true',
                        help='do not use FACS as a fallback for missing tissues')
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--output', default='data/mouse_atlas_cell_types.parquet')
    parser.add_argument('--age', default=DEFAULT_AGE, help='Tabula Muris Senis age (default: 3m)')
    args = parser.parse_args()

    input_file = args.input or DEFAULT_INPUT
    urls = utils.read_config(args.config, field='urls')
    if args.input is None:
        input_file = download_atlas(urls['tabula_muris_senis_droplet_url'], DEFAULT_INPUT)
    droplet = aggregate_expression(input_file=input_file, config_file=args.config, age=args.age)
    facs = None
    if not args.droplet_only:
        facs_file = args.facs_input or DEFAULT_FACS_INPUT
        if args.facs_input is None:
            facs_file = download_atlas(urls['tabula_muris_senis_facs_url'], DEFAULT_FACS_INPUT)
        facs = aggregate_expression(input_file=facs_file, config_file=args.config, age=args.age)
    data = combine_atlases(droplet, facs)
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    utils.save_to_parquet(data, args.output)
    print(f'Wrote {len(data):,} mouse cell-type expression rows to {args.output}')
