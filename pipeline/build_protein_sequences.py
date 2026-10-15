'''
Builds the protein sequence table the app uses for the FASTA downloads of the predicted
host-parasite PPI networks.
'''
import argparse
import io
import os
import sys

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils


def get_predicted_proteins(data_dir):
    '''Collects the proteins the app can show, grouped by species.'''
    predictions = utils.read_parquet_file(input_file=os.path.join(data_dir, 'predictions.parquet'))
    proteins = set(predictions['source']).union(predictions['target'])

    per_species = {}
    for protein in proteins:
        # STRING ids are <taxid>.<protein>, which is also how the sequence files are named
        per_species.setdefault(protein.split('.')[0], set()).add(protein)

    return per_species


def get_species_sequences(sequences_url, taxid, proteins):
    '''
    Parses the STRING protein.sequences file of one species and keeps the sequences of the
    requested proteins.
    '''
    filename = utils.download_file(url=sequences_url.replace('TAXID', str(taxid)),
                                   data_dir=os.path.join('data/downloads/species', str(taxid)))
    with utils.read_gzipped_file(filename) as handle:
        return [(record.id, str(record.seq))
                for record in SeqIO.parse(io.TextIOWrapper(handle), 'fasta')
                if record.id in proteins]


def build_sequences(data_dir, config_file, output_file):
    '''Writes the protein id --> sequence table of every protein in the predictions.'''
    sequences_url = utils.read_config(filepath=config_file, field='urls')['string_sequences_url']

    per_species = get_predicted_proteins(data_dir)
    sequences = []
    for taxid in sorted(per_species, key=int):
        species_sequences = get_species_sequences(sequences_url, taxid, per_species[taxid])
        print(f'    taxid {taxid}: {len(species_sequences)}/{len(per_species[taxid])} sequences found')
        sequences.extend(species_sequences)

    df = pd.DataFrame(sequences, columns=['protein', 'sequence'])
    utils.save_to_parquet(df, output_file)
    print(f'Wrote {len(df)} protein sequences to {output_file}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--output', default=None)
    args = parser.parse_args()

    output = args.output or os.path.join(args.data_dir, 'protein_sequences.parquet')
    build_sequences(args.data_dir, args.config, output)
