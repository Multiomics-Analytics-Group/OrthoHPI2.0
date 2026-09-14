'''
Builds the protein annotation table the app uses for the node tooltips of the predicted
host-parasite PPI networks.
'''
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils


# STRING writes "<name>; <function...>" and only the name belongs in a tooltip
MAX_DESCRIPTION_LENGTH = 120


def get_predicted_proteins(data_dir):
    '''Collects the proteins the app can show, grouped by species.'''
    predictions = utils.read_parquet_file(input_file=os.path.join(data_dir, 'predictions.parquet'))
    proteins = set(predictions['source']).union(predictions['target'])

    per_species = {}
    for protein in proteins:
        # STRING ids are <taxid>.<protein>, which is also how the info files are named
        per_species.setdefault(protein.split('.')[0], set()).add(protein)

    return per_species


def get_species_annotations(string_url, taxid, proteins):
    '''
    Parses the STRING protein.info file of one species and keeps the description of the
    requested proteins.
    '''
    filename = utils.download_file(url=string_url.replace('TAXID', str(taxid)),
                                   data_dir=os.path.join('data/downloads/species', str(taxid)))
    annotations = []
    with utils.read_gzipped_file(filename) as handle:
        next(handle, None)  # skip header
        for line in handle:
            fields = line.decode('utf-8').rstrip().split('\t')
            if len(fields) < 4 or fields[0] not in proteins:
                continue
            annotations.append((fields[0], shorten_annotation(fields[3])))

    return annotations


def shorten_annotation(annotation):
    '''
    Keeps the descriptive name at the head of a STRING annotation and drops the
    functional description that follows it, so the tooltip stays readable.
    '''
    description = annotation.split(';')[0].strip()
    if len(description) > MAX_DESCRIPTION_LENGTH:
        description = description[:MAX_DESCRIPTION_LENGTH].rstrip() + '...'

    return description


def build_annotations(data_dir, config_file, output_file):
    '''Writes the protein id --> description table of every protein in the predictions.'''
    string_url = utils.read_config(filepath=config_file, field='urls')['string_protein_url']

    per_species = get_predicted_proteins(data_dir)
    annotations = []
    for taxid in sorted(per_species, key=int):
        species_annotations = get_species_annotations(string_url, taxid, per_species[taxid])
        print(f'    taxid {taxid}: {len(species_annotations)}/{len(per_species[taxid])} proteins annotated')
        annotations.extend(species_annotations)

    df = pd.DataFrame(annotations, columns=['protein', 'description'])
    utils.save_to_parquet(df, output_file)
    print(f'Wrote {len(df)} protein descriptions to {output_file}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--output', default=None)
    args = parser.parse_args()

    output = args.output or os.path.join(args.data_dir, 'protein_annotations.parquet')
    build_annotations(args.data_dir, args.config, output)
