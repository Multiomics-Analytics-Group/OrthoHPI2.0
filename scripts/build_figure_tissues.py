'''
Annotate the host proteins of the predictions with the organs drawn on the TISSUES body
figures (images/tissues/tissues_<species>.svg).
'''
import argparse
import csv
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils
# the same host cutoff resolution the pipeline's tissue filter uses
from pipeline.filters import tissue_cutoff

# BTO code -> normalized organ name; source-label aliases are resolved in app/body_figure.py
FIGURE_ORGANS = {
    'BTO:0000047': 'adrenal gland',
    'BTO:0000089': 'blood',
    'BTO:0000140': 'bone',
    'BTO:0000141': 'bone marrow',
    'BTO:0000439': 'eye',
    'BTO:0000493': 'gall bladder',
    'BTO:0000562': 'heart',
    'BTO:0000648': 'intestine',
    'BTO:0000671': 'kidney',
    'BTO:0000759': 'liver',
    'BTO:0000763': 'lung',
    'BTO:0000784': 'lymph nodes',
    'BTO:0000887': 'muscle',
    'BTO:0000988': 'pancreas',
    'BTO:0001202': 'saliva',
    'BTO:0001253': 'skin',
    'BTO:0001281': 'spleen',
    'BTO:0001307': 'stomach',
    'BTO:0001379': 'thyroid gland',
    'BTO:0001418': 'urinary bladder',
    'BTO:0001419': 'urinary bladder',
    'BTO:0001484': 'nervous system',
}

# the pipeline's default TISSUES cutoff; --cutoff overrides it
DEFAULT_CUTOFF = 2.5

# columns of the Jensen Lab experiments file: protein, name, BTO code, label, source, source
# score, confidence
BTO_COLUMN = 2
SCORE_COLUMN = 6


def experiments_tissue_url(host):
    '''Return a host's experiments-channel URL from the shared pipeline configuration.'''
    url = host.get('tissues_url')
    if url is None:
        return None, None

    species = os.path.basename(url).split('_')[0]

    return species, url


def read_organ_annotations(annotation_file, taxid, valid_proteins, cutoff):
    '''
    Scan a Jensen Lab experiments tissue file for every organ the body figure draws that
    a protein passes the cutoff in.
    '''
    annotated = {}
    with open(annotation_file, 'r') as f:
        for data in csv.reader(f, delimiter='\t'):
            if len(data) <= SCORE_COLUMN:
                continue
            organ = FIGURE_ORGANS.get(data[BTO_COLUMN])
            if organ is None:
                continue
            protein = f'{taxid}.{data[0]}'
            if protein in valid_proteins:
                score = float(data[SCORE_COLUMN])
                if score >= cutoff:
                    # the same organ can be reached through more than one BTO code
                    annotated.setdefault(protein, set()).add(organ)

    return [(protein, organ) for protein, organs in annotated.items() for organ in organs]


def build(config_file, data_dir, cutoff):
    '''
    Annotates the host proteins of the predictions with the organs of the body figures.
    '''
    config = utils.read_config(filepath=config_file)
    predictions = utils.read_parquet_file(input_file=os.path.join(data_dir, 'predictions.parquet'))

    rows = []
    for taxid, host in config['hosts'].items():
        species, url = experiments_tissue_url(host)
        if species is None:
            print(f'  {taxid}: no tissues_url, skipped')
            continue

        # the annotation files cover the whole proteome
        valid_proteins = set(predictions.loc[predictions['taxid2'] == str(taxid), 'target'])
        if not valid_proteins:
            print(f'  {taxid} ({species}): no predictions, skipped')
            continue

        # the host's own cutoff, which decided the predictions being shaded
        host_cutoff = tissue_cutoff(host, cutoff)
        filename = utils.download_file(url=url, data_dir=os.path.join(data_dir, 'downloads'))
        host_rows = read_organ_annotations(filename, taxid, valid_proteins, host_cutoff)
        rows.extend(host_rows)
        annotated = len({protein for protein, _ in host_rows})
        print(f'  {taxid} ({species}): {annotated}/{len(valid_proteins)} host proteins '
              f'annotated at cutoff {host_cutoff}, {len(host_rows)} protein-organ pairs')

    df = pd.DataFrame(rows, columns=['Gene', 'Organ']).drop_duplicates()
    output_file = os.path.join(data_dir, 'figure_tissues.parquet')
    utils.save_to_parquet(df, output_file)
    print(f'wrote {output_file}: {len(df)} rows, {df["Organ"].nunique()} organs')

    return df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='config.yml', help='path to the configuration file')
    parser.add_argument('--data-dir', default='data', help='directory holding predictions.parquet')
    parser.add_argument('--cutoff', type=float, default=DEFAULT_CUTOFF,
                        help='minimum TISSUES confidence score for a host whose '
                             'configuration entry sets no tissue_cutoff')
    args = parser.parse_args()

    build(config_file=args.config, data_dir=args.data_dir, cutoff=args.cutoff)


if __name__ == '__main__':
    main()
