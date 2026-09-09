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


def get_parasite_niches(config_file, known_niches, default_niche):
    """
    Where each parasite of the configuration sits relative to the host cell, which decides
    which host proteins it is in a position to reach.

    :param str config_file: path to the configuration file
    :param known_niches: the niche names the caller has cut-offs for
    :param str default_niche: niche used for a parasite the config records none for
    :return: {parasite taxid (int): niche}
    """
    parasites = utils.read_config(filepath=config_file, field='parasites')
    niches = {}
    for taxid, parasite in parasites.items():
        niche = str(parasite.get('niche', '')).strip().lower()
        if niche not in known_niches:
            print(f"  WARNING: {parasite.get('label', taxid)} has no usable niche "
                  f"({parasite.get('niche')!r}); filtering it as {default_niche}")
            niche = default_niche
        niches[int(taxid)] = niche

    return niches


def host_niches(config_file, taxid, niches):
    """
    The niches that reach into one host: those of the parasites that infect it. A parasite
    with no `hosts` list is taken to infect every host, which is how get_links reads it.

    :param str config_file: path to the configuration file
    :param taxid: host taxid
    :param dict niches: {parasite taxid: niche}, from get_parasite_niches
    :return: set of niche names
    """
    parasites = utils.read_config(filepath=config_file, field='parasites')
    reaching = set()
    for parasite_taxid, parasite in parasites.items():
        hosts = parasite.get('hosts')
        if hosts is None or int(taxid) in hosts:
            reaching.add(niches[int(parasite_taxid)])

    return reaching


def apply_deeploc_filter(config_file, valid_proteins, deeploc_dir, niche_cutoffs,
                         default_niche):
    """
    Filter host proteins to the ones a parasite can reach, using DeepLoc 2 (Accurate)
    predictions, replacing the COMPARTMENTS plasma-membrane filter.

    Which localisation classes count depends on the parasite's niche, so the filter has
    two halves. valid_proteins is a pool per host and not per parasite, so it is narrowed
    here to the union of the classes any parasite infecting that host reaches -- the
    surface for a host met by extracellular parasites alone, the surface plus the cytosol
    and nucleus for one that also carries an intracellular parasite. The returned sets are
    the per-niche halves of that pool, which homology.get_links applies parasite by
    parasite so that an extracellular parasite is not handed a cytosolic host protein.

    A host protein is kept for a niche if any of the niche's classes scores above its
    cut-off. Strictly greater than, which is how DeepLoc itself calls a class
    (DeepLoc2/utils.py convert_label2string), and the same comparison the parasite
    secretome filter makes in deeploc/build_secretome_fastas.py.

    The probability is read rather than the Localizations column DeepLoc writes beside it:
    that column names a class even for a protein that crosses no threshold at all, falling
    back to whichever came closest, which is not evidence of anything.

    DeepLoc Protein_IDs are already STRING ids (taxid.<protein>), so they match the
    valid_proteins keys directly.

    :param str config_file: path to the configuration file
    :param dict valid_proteins: {taxid: {protein_id: name}}; each host is filtered in place
    :param str deeploc_dir: directory of DeepLoc results (<deeploc_dir>/<taxid>/results_*.csv)
    :param dict niche_cutoffs: {niche: {DeepLoc class: probability cut-off}}
    :param str default_niche: niche for a parasite the config records none for
    :return: {niche: set of host protein ids it reaches}, over every host at once
    """
    hosts = utils.read_config(filepath=config_file, field='hosts')
    niches = get_parasite_niches(config_file, set(niche_cutoffs), default_niche)
    reachable = {niche: set() for niche in niche_cutoffs}
    for taxid in hosts:
        reaching = host_niches(config_file, taxid, niches)
        matches = sorted(glob.glob(os.path.join(deeploc_dir, str(taxid), 'results_*.csv')))
        if not matches:
            print(f"  WARNING: no DeepLoc results for host {taxid} in {deeploc_dir}; "
                  "keeping its proteins unfiltered")
            for niche in reaching:
                reachable[niche].update(valid_proteins[taxid])
            continue

        df = pd.read_csv(matches[-1])
        per_niche = {}
        for niche in reaching:
            keep = pd.Series(False, index=df.index)
            for localisation, cutoff in niche_cutoffs[niche].items():
                keep = keep | (df[localisation] > cutoff)
            per_niche[niche] = set(df.loc[keep, 'Protein_ID'])

        kept = set().union(*per_niche.values()) if per_niche else set()
        valid_proteins[taxid] = {p: n for p, n in valid_proteins[taxid].items()
                                 if p in kept}
        for niche, proteins in per_niche.items():
            reachable[niche].update(proteins & set(valid_proteins[taxid]))
        print(f"    taxid {taxid}: {len(valid_proteins[taxid])} proteins kept ("
              + ', '.join(f"{niche} {len(reachable[niche] & set(valid_proteins[taxid]))}"
                          for niche in sorted(per_niche)) + ")")

    return reachable


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
