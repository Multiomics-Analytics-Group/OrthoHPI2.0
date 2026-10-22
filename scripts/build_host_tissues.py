'''
Write host_tissues.parquet: every host protein the TISSUES filter passes, with the
tissues it is annotated to, before the DeepLoc filter narrows them. The hosts-of-a-parasite
page reads it to tell a group that is not expressed in an infected tissue from one that
DeepLoc puts out of reach; tissues_cell_types.parquet holds only proteins past both.
'''
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from pipeline import filters, main


def build(config_file, data_dir):
    print(f"Getting proteins for {config_file}...")
    proteins = main.get_proteins(config_file)

    print("Applying the tissue filter...")
    tissues = filters.apply_tissue_filter(config_file=config_file, valid_proteins=proteins,
                                          cutoff=main.TISSUE_CUTOFF)
    table = pd.DataFrame([(protein, tissue) for protein, labels in tissues.items()
                          for tissue in labels], columns=['protein', 'tissue'])
    print(f"  {table['protein'].nunique()} host proteins in an annotated tissue")

    output_file = os.path.join(data_dir, 'host_tissues.parquet')
    utils.save_to_parquet(table, output_file)
    print(f"Wrote {output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml', help='configuration file to read')
    parser.add_argument('--data-dir', default='data', help='directory to write into')
    args = parser.parse_args()

    build(config_file=args.config, data_dir=args.data_dir)
