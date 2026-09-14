'''
Write proteome_sizes.parquet: how many proteins STRING holds for every species of the
config, before any filter. The home page draws it against the eligible pool.
'''
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pipeline import main


def build(config_file, data_dir):
    print(f"Getting proteins for {config_file}...")
    proteins = main.get_proteins(config_file)
    sizes = pd.DataFrame({'taxid': [str(taxid) for taxid in proteins],
                          'proteins': [len(names) for names in proteins.values()]})

    output_file = os.path.join(data_dir, 'proteome_sizes.parquet')
    sizes.to_parquet(output_file, index=False)
    print(f"Wrote {output_file}: {len(sizes)} species, "
          f"{sizes['proteins'].sum():,} proteins")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml', help='configuration file to read')
    parser.add_argument('--data-dir', default='data', help='directory to write into')
    args = parser.parse_args()

    build(config_file=args.config, data_dir=args.data_dir)
