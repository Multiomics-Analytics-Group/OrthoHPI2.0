'''
Builds the DeepLoc localisation table the app uses to say where a predicted interactor sits
in the cell.
'''
import argparse
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils


# the DeepLoc results directory the pipeline filters on, relative to the data directory
DEEPLOC_ACCURATE_DIR = os.path.join('deeploc', 'output_accurate', 'deeploc_output_accurate')
# the four classes the host filter reads (pipeline/main.py DEEPLOC_NICHE_CLASSES) and the
# three text columns
COLUMNS = {'Protein_ID': 'protein', 'Localizations': 'localizations', 'Signals': 'signals',
           'Membrane types': 'membrane_types', 'Extracellular': 'extracellular',
           'Cell membrane': 'cell_membrane', 'Cytoplasm': 'cytoplasm',
           'Nucleus': 'nucleus'}


def get_predicted_proteins(data_dir):
    '''Collects the proteins the app can show, grouped by species.'''
    predictions = utils.read_parquet_file(input_file=os.path.join(data_dir, 'predictions.parquet'))
    proteins = set(predictions['source']).union(predictions['target'])

    per_species = {}
    for protein in proteins:
        # STRING ids are <taxid>.<protein>, which is also how the results are named
        per_species.setdefault(protein.split('.')[0], set()).add(protein)

    return per_species


def get_species_localisations(deeploc_dir, taxid, proteins):
    '''Reads the DeepLoc results of one species and keeps the requested proteins.'''
    matches = sorted(glob.glob(os.path.join(deeploc_dir, str(taxid), 'results_*.csv')))
    if not matches:
        return None

    df = pd.read_csv(matches[-1], usecols=list(COLUMNS))
    df = df[df['Protein_ID'].isin(proteins)].rename(columns=COLUMNS)

    return df[list(COLUMNS.values())]


def build_localisations(data_dir, deeploc_dir, output_file):
    '''
    Writes the localisation table of every protein in the predictions DeepLoc was run
    on.
    '''
    per_species = get_predicted_proteins(data_dir)
    localisations = []
    for taxid in sorted(per_species, key=int):
        df = get_species_localisations(deeploc_dir, taxid, per_species[taxid])
        if df is None:
            print(f'    taxid {taxid}: no DeepLoc results in {deeploc_dir}, skipped')
            continue
        print(f'    taxid {taxid}: {len(df)}/{len(per_species[taxid])} proteins localised')
        localisations.append(df)

    df = pd.concat(localisations, ignore_index=True)
    # an empty string reads better in a tooltip than NaN
    for column in ['localizations', 'signals', 'membrane_types']:
        df[column] = df[column].fillna('')
    utils.save_to_parquet(df, output_file)
    print(f'Wrote {len(df)} protein localisations to {output_file}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--deeploc-dir', default=None)
    parser.add_argument('--output', default=None)
    args = parser.parse_args()

    deeploc = args.deeploc_dir or os.path.join(args.data_dir, DEEPLOC_ACCURATE_DIR)
    output = args.output or os.path.join(args.data_dir, 'deeploc_localisations.parquet')
    build_localisations(args.data_dir, deeploc, output)
