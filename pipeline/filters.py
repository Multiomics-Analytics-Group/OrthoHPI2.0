import glob
import os

import pandas as pd

import utils


def tissue_cutoff(host, default):
    """Return a host-specific TISSUES cutoff, or the pipeline default when absent."""
    cutoff = host.get('tissue_cutoff', default)
    if not isinstance(cutoff, (int, float)):
        raise TypeError(f"tissue_cutoff for {host.get('label', 'host')} must be numeric")
    return float(cutoff)


def apply_tissue_filter(config_file, valid_proteins, cutoff):
    hosts = utils.read_config(filepath=config_file, field='hosts')
    tissue_mapping = utils.read_config(filepath=config_file, field='tissues')
    tissues = {}
    for taxid in hosts:
        proteins = valid_proteins[taxid]
        if 'tissues_url' in hosts[taxid]:
            url = hosts[taxid]['tissues_url']
            filename = utils.download_file(url=url, data_dir='data/downloads')
            host_tissues, proteins = get_tissues(
                config_file, filename, proteins, tissue_cutoff(hosts[taxid], cutoff),
                tissue_mapping, taxid)
            tissues.update(host_tissues)

        valid_proteins[taxid] = proteins

    return tissues


def _filter_by_annotation(annotation_file, valid_proteins, cutoff, taxid, valid_values, score_col, transform):
    """
    Scan a jensenlab TSV and keep proteins whose annotation (column 2) is in valid_values
    with a confidence score (score_col) of at least cutoff.

    :param str annotation_file: path to the jensenlab tissue/compartment TSV
    :param dict valid_proteins: {protein_id: name} candidates for this host
    :param float cutoff: minimum confidence score accepted
    :param taxid: host taxid, used to prefix the file's protein ids into STRING ids
    :param set valid_values: annotation values (BTO tissue / GO compartment) to keep
    :param int score_col: column index holding the confidence score
    :param callable transform: maps a raw annotation value to how it is stored
    :return: (annotations, kept) where annotations is {protein: [transformed value, ...]}
             and kept is the {protein: name} subset that passed the filter
    """
    annotations = {}
    kept = {}
    with open(annotation_file, 'r') as f:
        next(f, None)  # skip header
        for line in f:
            data = line.rstrip().split('\t')
            protein = f"{taxid}." + data[0]
            value = data[2]
            score = float(data[score_col])
            if protein in valid_proteins and score >= cutoff and value in valid_values:
                annotations.setdefault(protein, []).append(transform(value))
                kept[protein] = valid_proteins[protein]

    return annotations, kept


def get_tissues(config_file, tissues_file, valid_proteins, cutoff, mapping, taxid):
    """
    Get protein tissue expression for tissues relevant in the lifecycle of the
    studied parasites.

    :param str config_file: path to the configuration file
    :param str tissues_file: path to file with tissue expression (tissues.jensenlab.org)
    :param dict valid_proteins: {protein_id: name} candidates for this host
    :param float cutoff: minimum confidence score accepted (tissues.jensenlab.org)
    :param dict mapping: BTO tissue code -> display name
    :param taxid: host taxid, used to prefix protein ids into STRING ids
    :return: (tissues, kept) — {protein: [tissue name, ...]} and the passing {protein: name}
    """
    parasites = utils.read_config(filepath=config_file, field='parasites')
    valid_tissues = set()
    for parasite in parasites.values():
        valid_tissues.update(parasite['tissues'])

    return _filter_by_annotation(tissues_file, valid_proteins, cutoff, taxid,
                                 valid_tissues, score_col=6, transform=lambda t: mapping[t])


def apply_deeploc_filter(config_file, valid_proteins, deeploc_dir,
                         extracellular_cutoff, membrane_cutoff):
    """
    Filter host proteins to surface-exposed ones using DeepLoc 2 (Accurate) predictions,
    replacing the COMPARTMENTS plasma-membrane filter.

    A host protein is kept if P(Cell membrane) >= membrane_cutoff, or (when
    extracellular_cutoff is not None) also if P(Extracellular) >= extracellular_cutoff.
    DeepLoc Protein_IDs are already STRING ids (taxid.<protein>), so they match the
    valid_proteins keys directly.

    :param str config_file: path to the configuration file
    :param dict valid_proteins: {taxid: {protein_id: name}}; each host is filtered in place
    :param str deeploc_dir: directory of DeepLoc results (<deeploc_dir>/<taxid>/results_*.csv)
    :param extracellular_cutoff: minimum P(Extracellular) to keep a protein, or None to
                                 keep only Cell membrane proteins
    :param float membrane_cutoff: minimum P(Cell membrane) to keep a protein
    :return: valid_proteins with non-surface host proteins removed
    """
    hosts = utils.read_config(filepath=config_file, field='hosts')
    for taxid in hosts:
        matches = sorted(glob.glob(os.path.join(deeploc_dir, str(taxid), 'results_*.csv')))
        if not matches:
            print(f"  WARNING: no DeepLoc results for host {taxid} in {deeploc_dir}; "
                  "keeping its proteins unfiltered")
            continue

        df = pd.read_csv(matches[-1])
        keep = df['Cell membrane'] >= membrane_cutoff
        if extracellular_cutoff is not None:
            keep = keep | (df['Extracellular'] >= extracellular_cutoff)
        surface_ids = set(df.loc[keep, 'Protein_ID'])
        valid_proteins[taxid] = {p: n for p, n in valid_proteins[taxid].items()
                                 if p in surface_ids}

    return valid_proteins


def get_secretome_predictions(config_file, secretome_dir, valid_proteins):
    """
    Filter out proteins that are not secreted or membrane from the list of parasite proteins

    :param str config_file: path to the configuration file
    :param str secretome_dir: path to the directory where the prediction files are
    :param dict valid_proteins: dictionary with annotations in valid proteins

    :return filtered_dict: dictionary with only secreted or membrane parasite proteins
    """
    parasites = utils.read_config(filepath=config_file, field='parasites')
    for parasite in parasites:
        filepath = os.path.join(secretome_dir, str(parasite)+'.fasta')
        sequences = utils.read_fasta(filepath)
        filter_out_ids = utils.filter_sequences(sequences, valid_proteins[parasite])
        for k in filter_out_ids:
            valid_proteins[parasite].pop(k, None)

    return valid_proteins
