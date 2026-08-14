"""
Builds the DeepLoc localisation table the app uses to say where a predicted
interactor sits in the cell.

DeepLoc is run outside the pipeline (see docs/deeploc.md) and the main pipeline
only ever reads its results as a filter: a host protein is kept if it is called
surface-exposed, and the probabilities behind that call are then thrown away.
This script keeps them for the proteins that survived into the predictions and
writes them to a small parquet the app reads alongside the predictions, so a
figure can show whether a host protein is a membrane protein or a secreted one
rather than only that it passed the filter.

Both sides are written: the host proteins the filter selected and the parasite
proteins, whose secretome filter is a different prediction (SignalP/TargetP) but
which DeepLoc was run on all the same.

Run it after the main pipeline, from the repo root:

    .venv/bin/python pipeline/build_deeploc_localisations.py

    --data-dir     directory holding predictions.parquet (default: data)
    --deeploc-dir  DeepLoc results, <dir>/<taxid>/results_*.csv
                   (default: <data-dir>/deeploc/output_accurate/deeploc_output_accurate)
    --output       parquet to write (default: <data-dir>/deeploc_localisations.parquet)
"""
import argparse
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils


# the DeepLoc results directory, relative to the data directory. The same one the
# pipeline filters on (pipeline/main.py: DEEPLOC_ACCURATE_DIR)
DEEPLOC_ACCURATE_DIR = os.path.join('deeploc', 'output_accurate', 'deeploc_output_accurate')
# columns kept out of the results file: the two surface classes the host filter is made
# of, and the three text columns naming what DeepLoc called the protein. The other eight
# compartment probabilities are of no use to a page about surface interactions
COLUMNS = {'Protein_ID': 'protein', 'Localizations': 'localizations', 'Signals': 'signals',
           'Membrane types': 'membrane_types', 'Extracellular': 'extracellular',
           'Cell membrane': 'cell_membrane'}


def get_predicted_proteins(data_dir):
    """
    Collects the proteins the app can show, grouped by species.

    :param str data_dir: directory holding predictions.parquet
    :return: {taxid as str: set of STRING protein ids}
    """
    predictions = utils.read_parquet_file(input_file=os.path.join(data_dir, 'predictions.parquet'))
    proteins = set(predictions['source']).union(predictions['target'])

    per_species = {}
    for protein in proteins:
        # STRING ids are <taxid>.<protein>, which is also how the results are named
        per_species.setdefault(protein.split('.')[0], set()).add(protein)

    return per_species


def get_species_localisations(deeploc_dir, taxid, proteins):
    """
    Reads the DeepLoc results of one species and keeps the requested proteins. The
    Protein_IDs are already STRING ids, so they match the predictions directly.

    Where a species was run more than once the newest results file is the one read,
    which is what the pipeline filter does with the same directory.

    :param str deeploc_dir: directory of DeepLoc results (<deeploc_dir>/<taxid>/results_*.csv)
    :param str taxid: taxonomic id of the species of interest
    :param set proteins: STRING protein ids to keep
    :return: dataframe of the kept proteins, or None if the species was never run
    """
    matches = sorted(glob.glob(os.path.join(deeploc_dir, str(taxid), 'results_*.csv')))
    if not matches:
        return None

    df = pd.read_csv(matches[-1], usecols=list(COLUMNS))
    df = df[df['Protein_ID'].isin(proteins)].rename(columns=COLUMNS)

    return df[list(COLUMNS.values())]


def build_localisations(data_dir, deeploc_dir, output_file):
    """
    Writes the localisation table of every protein in the predictions DeepLoc was run on.

    :param str data_dir: directory holding predictions.parquet
    :param str deeploc_dir: directory of DeepLoc results
    :param str output_file: parquet path to write
    """
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
    # the text columns are empty for a protein DeepLoc called no signal or membrane type,
    # and an empty string reads better in a tooltip than the NaN pandas gives it
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
