'''
Write filter_stages.parquet: per host-parasite pair, how many proteins of either side
survive each filter of the pipeline, counted as scripts/build_filter_funnel.py counts
them. The home page draws it as a funnel, and counts the last stage, the proteins in a
predicted interaction, itself at the confidence the page is set to.
'''
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from scripts.build_filter_funnel import count_pairs, load


def build(config_file, data_dir):
    config = utils.read_config(filepath=config_file)
    pairs = count_pairs(config, load(config_file, data_dir))
    # count_pairs names a host by its common name; the page reads it by its label
    labels = {h['label'].split('(')[1].rstrip(')'): h['label'] for h in config['hosts'].values()}
    pairs['host'] = pairs['host'].map(labels)
    pairs = pairs.drop(columns=['host_interactor', 'parasite_interactor', 'edges'])

    output_file = os.path.join(data_dir, 'filter_stages.parquet')
    pairs.to_parquet(output_file, index=False)
    print(f"Wrote {output_file}: {len(pairs)} pairs")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml', help='configuration file to read')
    parser.add_argument('--data-dir', default='data', help='directory to write into')
    args = parser.parse_args()

    build(config_file=args.config, data_dir=args.data_dir)
