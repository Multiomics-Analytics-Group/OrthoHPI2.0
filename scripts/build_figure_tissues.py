"""
Annotate the host proteins of the predictions with the organs drawn on the TISSUES
body figures (images/tissues/tissues_<species>.svg).

This is deliberately separate from the tissue annotation the pipeline already builds.
`config.yml` picks 33 fine-grained BTO terms chosen for parasite lifecycles, which is
what decides whether a prediction is relevant (pipeline/filters.py, filter_tissues in
app/web_utils.py). The body figure draws a different, coarser vocabulary: 21 organs,
which TISSUES also annotates and which the download already carries as its own rows --
BTO propagates a protein's annotation up to the parent organ, so a brain protein comes
with `nervous system` too. Reading those 21 codes straight out of the experiments file
used by the pipeline lets the figure be shaded without mapping the 33 lifecycle names
onto it by hand. This keeps the figure's annotations on the same direct-expression
evidence standard as predictions; organs without experimental rows remain unshaded.

Writes <data_dir>/figure_tissues.parquet with one row per (host protein, organ).
"""
import argparse
import csv
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils

# BTO code -> the normalized organ name. The source figures call the urinary-bladder
# region `urine` and differ on the thyroid label; these source-label aliases are resolved
# when the figure is drawn (app/body_figure.py).
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

# Minimum confidence score for figure annotations. --cutoff overrides it.
DEFAULT_CUTOFF = 2.5

# Organs kept per host protein, best score first.
DEFAULT_TOP_ORGANS = 3

# Columns of the Jensen Lab experiments file: protein, name, BTO code, label, source,
# source-specific score, confidence score.
BTO_COLUMN = 2
SCORE_COLUMN = 6


def experiments_tissue_url(host):
    """
    Return a host's experiments-channel URL from the shared pipeline configuration.

    :param dict host: one entry of the `hosts` section of the configuration
    :return: (species name as jensenlab spells it, download url), or (None, None) when
             the host has no tissue annotation
    """
    url = host.get('tissues_url')
    if url is None:
        return None, None

    species = os.path.basename(url).split('_')[0]

    return species, url


def read_organ_annotations(annotation_file, taxid, valid_proteins, cutoff, top_organs):
    """
    Scan a Jensen Lab experiments tissue file for the organs the body figure draws, keeping
    the best-scoring ones of each protein.

    The files run to hundreds of megabytes, so they are streamed rather than read into a
    dataframe, and only the proteins the predictions actually contain are kept.

    :param str annotation_file: path to the <species>_tissue_experiments_full.tsv
    :param taxid: host taxid, used to prefix the file's protein ids into STRING ids
    :param set valid_proteins: STRING ids to keep
    :param float cutoff: minimum confidence score accepted
    :param int top_organs: organs kept per protein, best score first
    :return: list of (protein, organ) pairs
    """
    scored = {}
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
                    # the same organ can be reached through more than one BTO code, and
                    # the pig figure labels one of them differently, so the best wins
                    previous = scored.setdefault(protein, {}).get(organ)
                    if previous is None or score > previous:
                        scored[protein][organ] = score

    rows = []
    for protein, organs in scored.items():
        best = sorted(organs, key=lambda organ: (-organs[organ], organ))[:top_organs]
        rows.extend((protein, organ) for organ in best)

    return rows


def build(config_file, data_dir, cutoff, top_organs):
    config = utils.read_config(filepath=config_file)
    predictions = utils.read_parquet_file(input_file=os.path.join(data_dir, 'predictions.parquet'))

    rows = []
    for taxid, host in config['hosts'].items():
        species, url = experiments_tissue_url(host)
        if species is None:
            print(f'  {taxid}: no tissues_url, skipped')
            continue

        # the host proteins the predictions target; the annotation files cover the whole
        # proteome, of which the network only ever shows a few hundred proteins
        valid_proteins = set(predictions.loc[predictions['taxid2'] == str(taxid), 'target'])
        if not valid_proteins:
            print(f'  {taxid} ({species}): no predictions, skipped')
            continue

        filename = utils.download_file(url=url, data_dir=os.path.join(data_dir, 'downloads'))
        host_rows = read_organ_annotations(filename, taxid, valid_proteins, cutoff, top_organs)
        rows.extend(host_rows)
        annotated = len({protein for protein, _ in host_rows})
        print(f'  {taxid} ({species}): {annotated}/{len(valid_proteins)} host proteins '
              f'annotated, {len(host_rows)} protein-organ pairs')

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
                        help='minimum TISSUES confidence score')
    parser.add_argument('--top-organs', type=int, default=DEFAULT_TOP_ORGANS,
                        help='organs kept per host protein, best score first')
    args = parser.parse_args()

    build(config_file=args.config, data_dir=args.data_dir, cutoff=args.cutoff,
          top_organs=args.top_organs)


if __name__ == '__main__':
    main()
