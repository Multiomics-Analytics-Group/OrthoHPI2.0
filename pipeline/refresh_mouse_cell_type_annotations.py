"""Add generated mouse atlas annotations to the existing app tissue artifact."""
import argparse
import os

import pandas as pd

import utils
from . import cell_type_annotations

MOUSE_TAXID = f'{cell_type_annotations.MOUSE_TAXID}.'


def refresh(data_dir):
    """Replace mouse cell-type rows in tissues_cell_types.parquet from the atlas artifact."""
    output_file = os.path.join(data_dir, 'tissues_cell_types.parquet')
    tissues = utils.read_parquet_file(output_file)
    required = {'Gene', 'Tissue'}
    missing = required.difference(tissues.columns)
    if missing:
        raise ValueError(f'{output_file} is missing required columns: {sorted(missing)}')

    is_mouse = tissues['Gene'].astype(str).str.startswith(MOUSE_TAXID)
    mouse_tissues = tissues.loc[is_mouse, ['Gene', 'Tissue']].drop_duplicates()
    mouse_data = cell_type_annotations.read_mouse_atlas(
        data_dir=data_dir, valid_proteins=mouse_tissues['Gene'])
    mouse = mouse_tissues.merge(mouse_data, on=['Gene', 'Tissue'], how='left')
    data = pd.concat([tissues.loc[~is_mouse], mouse], ignore_index=True, sort=False)
    utils.save_to_parquet(data, output_file)
    print(f'Updated mouse cell-type annotations in {output_file}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', default='data')
    args = parser.parse_args()
    refresh(data_dir=args.data_dir)
