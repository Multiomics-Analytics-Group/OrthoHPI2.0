"""
Convert cell-type means extracted from the pig atlas into the app's Parquet format.

The raw ``pig_atlas_20221014.rds`` Seurat object is too large to load as part of
the prediction pipeline. First aggregate its normalized sparse expression matrix:

    Rscript scripts/extract_pig_atlas_expression.R data/pig_atlas_20221014.rds \
        /tmp/pig_atlas_expression.csv

Then run this script from the repository root:

    .venv/bin/python -m pipeline.build_pig_atlas_cell_types \
        --input /tmp/pig_atlas_expression.csv

The result is ``data/pig_atlas_cell_types.parquet``. On the next pipeline run,
``pipeline.main`` merges it with human HPA annotation when writing
``tissues_cell_types.parquet``.
"""
import argparse
import os

import pandas as pd

import utils

PIG_TAXID = '9823'
REQUIRED_COLUMNS = {'Gene', 'Tissue', 'Cell type', 'nTPM'}
TISSUE_MAPPING = {
    'Intestine': 'intestine',
    'Neonatal Ileum_D0': 'intestine',
    'Neonatal Ileum_D01': 'intestine',
    'Neonatal Ileum_D03': 'intestine',
    'Neonatal Ileum_D07': 'intestine',
    'Neonatal Ileum_D14': 'intestine',
    'Neonatal Ileum_D21': 'intestine',
    'PBMC': 'blood',
}


def map_expression(input_file, config_file):
    """Map pig-atlas genes and tissues to STRING IDs and OrthoHPI tissue names."""
    data = pd.read_csv(input_file)
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        raise ValueError(f'{input_file} is missing required columns: {sorted(missing)}')

    aliases = utils.parse_string_aliases(config_file=config_file,
                                         sources=['Ensembl_gene'], taxid=PIG_TAXID)
    valid_tissues = {name.lower() for name in
                     utils.read_config(filepath=config_file, field='tissues').values()}

    data = data[list(REQUIRED_COLUMNS)].copy()
    data['Gene'] = data['Gene'].map(aliases)
    data['Tissue'] = data['Tissue'].replace(TISSUE_MAPPING).str.lower()
    data['nTPM'] = pd.to_numeric(data['nTPM'], errors='raise')
    data = data[data['Gene'].notna() & data['Tissue'].isin(valid_tissues) & (data['nTPM'] > 0)]
    data = (data.groupby(['Gene', 'Tissue', 'Cell type'], as_index=False, observed=True)['nTPM']
            .mean())
    return data


def build(input_file, config_file, output_file):
    """Write mapped pig cell-type expression in the shared annotation schema."""
    data = map_expression(input_file=input_file, config_file=config_file)
    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    utils.save_to_parquet(data, output_file)
    print(f'Wrote {len(data):,} pig cell-type expression rows to {output_file}')
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', required=True,
                        help='CSV produced by scripts/extract_pig_atlas_expression.R')
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--output', default='data/pig_atlas_cell_types.parquet')
    args = parser.parse_args()
    build(input_file=args.input, config_file=args.config, output_file=args.output)
