"""
Builds the protein annotation table the app uses for the node tooltips of the
predicted host-parasite PPI networks.

The STRING protein.info files already downloaded per species carry a descriptive
name for every protein ("Beta-2-microglobulin; Component of the class I major
histocompatibility complex..."), which the predictions table does not keep -- it
only keeps the short preferred name shown as the node label. This script pulls
those descriptions for the proteins that appear in the predictions and writes
them to a small parquet the app reads alongside the predictions.

Run it after the main pipeline, from the repo root:

    .venv/bin/python pipeline/build_protein_annotations.py

    --data-dir    directory holding predictions.parquet (default: data)
    --config      configuration file with the STRING urls (default: config.yml)
    --output      parquet to write (default: <data-dir>/protein_annotations.parquet)
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils


# the description is the whole annotation up to the first sentence break; STRING
# writes "<name>; <function...>" and only the name belongs in a tooltip
MAX_DESCRIPTION_LENGTH = 120


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
        # STRING ids are <taxid>.<protein>, which is also how the info files are named
        per_species.setdefault(protein.split('.')[0], set()).add(protein)

    return per_species


def get_species_annotations(string_url, taxid, proteins):
    """
    Parses the STRING protein.info file of one species and keeps the description of
    the requested proteins. The file is downloaded if it is not in data/downloads yet.

    :param str string_url: url template to the STRING protein.info file (contains TAXID)
    :param str taxid: taxonomic id of the species of interest
    :param set proteins: STRING protein ids to keep
    :return: list of (protein id, description) tuples
    """
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
    """
    Keeps the descriptive name at the head of a STRING annotation and drops the
    functional description that follows it, so the tooltip stays readable.

    :param str annotation: annotation field of the STRING protein.info file
    :return: the descriptive name of the protein
    """
    description = annotation.split(';')[0].strip()
    if len(description) > MAX_DESCRIPTION_LENGTH:
        description = description[:MAX_DESCRIPTION_LENGTH].rstrip() + '...'

    return description


def build_annotations(data_dir, config_file, output_file):
    """
    Writes the protein id --> description table of every protein in the predictions.

    :param str data_dir: directory holding predictions.parquet
    :param str config_file: path to the configuration file
    :param str output_file: parquet path to write
    """
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
