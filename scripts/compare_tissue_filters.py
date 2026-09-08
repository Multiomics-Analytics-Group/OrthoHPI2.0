"""
Compare how much the jensenlab TISSUES filter and an HPA-expression filter each
narrow the host protein pool, and how much they agree.

Standalone analysis -- it reuses the pipeline's own functions but writes nothing
the pipeline reads.

The jensenlab channel can be swapped: `experiments` is what the pipeline uses,
`integrated` folds in text-mining and curated knowledge on top of it. The two files
have the same layout apart from where the confidence score sits.

Usage: .venv/bin/python scripts/compare_tissue_filters.py [experiments|integrated]
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from pipeline import cell_type_annotations, filters, main

CONFIG = 'config.yml'
HPA_NTPM_CUTOFFS = [0.0, 1.0]

# Column holding the confidence score in each jensenlab channel's TSV.
SCORE_COL = {'experiments': 6, 'integrated': 4}
SCORE_SWEEP = [0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]


def host_pools(config_file):
    """{taxid: {protein: name}} straight out of STRING, hosts only."""
    hosts = utils.read_config(filepath=config_file, field='hosts')
    urls = utils.read_config(filepath=config_file, field='urls')
    return {taxid: main.get_species_proteins(urls['string_protein_url'], taxid)
            for taxid in hosts}


def tissue_file(config_file, taxid, channel):
    """Download the host's jensenlab TSV for the requested channel."""
    hosts = utils.read_config(filepath=config_file, field='hosts')
    url = hosts[taxid]['tissues_url'].replace('experiments', channel)
    return utils.download_file(url=url, data_dir='data/downloads')


def parasite_tissues(config_file):
    """The BTO tissues some parasite of the config infects -- what the filter scans for."""
    parasites = utils.read_config(filepath=config_file, field='parasites')
    valid_tissues = set()
    for parasite in parasites.values():
        valid_tissues.update(parasite['tissues'])
    return valid_tissues


def best_scores(config_file, pool, taxid, channel):
    """{protein: best confidence score} over the tissues the config's parasites infect."""
    valid_tissues = parasite_tissues(config_file)
    score_col = SCORE_COL[channel]
    best = {}
    with open(tissue_file(config_file, taxid, channel)) as f:
        for line in f:
            data = line.rstrip().split('\t')
            protein = f"{taxid}." + data[0]
            if protein in pool[taxid] and data[2] in valid_tissues:
                score = float(data[score_col])
                if score > best.get(protein, -1):
                    best[protein] = score
    return best


def tissues_pass(config_file, pool, cutoff, channel):
    """Proteins kept by the jensenlab TISSUES filter, per host."""
    hosts = utils.read_config(filepath=config_file, field='hosts')
    mapping = utils.read_config(filepath=config_file, field='tissues')
    kept = {}
    for taxid in hosts:
        _, proteins = filters._filter_by_annotation(
            tissue_file(config_file, taxid, channel), dict(pool[taxid]), cutoff, taxid,
            parasite_tissues(config_file), SCORE_COL[channel],
            transform=lambda t: mapping[t])
        kept[taxid] = set(proteins)
    return kept


def deeploc_pass(config_file, pool):
    """Proteins kept by the DeepLoc surface filter, per host."""
    proteins = {taxid: dict(p) for taxid, p in pool.items()}
    filters.apply_deeploc_filter(
        config_file=config_file, valid_proteins=proteins,
        deeploc_dir=os.path.join('data', main.DEEPLOC_ACCURATE_DIR),
        extracellular_cutoff=main.DEEPLOC_EXTRACELLULAR_CUTOFF,
        membrane_cutoff=main.DEEPLOC_MEMBRANE_CUTOFF)
    return {taxid: set(p) for taxid, p in proteins.items()}


def hpa_pass(config_file, cutoffs):
    """Human proteins with HPA single-cell expression in a config tissue, per nTPM cutoff."""
    data = cell_type_annotations.read_hpa(config_file=config_file)
    mapped = cell_type_annotations.map_hpa_data(config_file=config_file, hpa_data=data)
    return {c: set(mapped.loc[mapped['nTPM'] > c, 'Gene'].dropna()) for c in cutoffs}


def hpa_any_tissue_universe(config_file):
    """Human STRING proteins HPA has a row for in ANY tissue (not just config tissues)."""
    data = cell_type_annotations.read_hpa(config_file=config_file)
    aliases = utils.parse_string_aliases(config_file, sources=['Ensembl_gene'])
    return set(data['Gene'].map(aliases).dropna())


def pct(n, d):
    return f"{100.0 * n / d:5.1f}%" if d else "    n/a"


def main_(channel='experiments'):
    pool = host_pools(CONFIG)
    hosts = utils.read_config(filepath=CONFIG, field='hosts')

    print(f"jensenlab channel = {channel}, TISSUES cutoff = {main.TISSUE_CUTOFF}\n")
    tis = tissues_pass(CONFIG, pool, main.TISSUE_CUTOFF, channel)
    dl = deeploc_pass(CONFIG, pool)

    print("=== Score sweep: proteins passing each cutoff, per host ===")
    print(f"{'host':<22}{'STRING':>9}" + ''.join(f"{'>=' + str(c):>10}" for c in SCORE_SWEEP))
    scores = {taxid: best_scores(CONFIG, pool, taxid, channel) for taxid in hosts}
    for taxid, meta in hosts.items():
        total = len(pool[taxid])
        counts = [sum(1 for v in scores[taxid].values() if v >= c) for c in SCORE_SWEEP]
        print(f"{meta['label']:<22}{total:>9}" + ''.join(f"{n:>10}" for n in counts))
        print(f"{'':<22}{'':>9}" + ''.join(f"{pct(n, total):>10}" for n in counts))
    print()

    print("=== Per host: how much each filter removes (from the full STRING pool) ===")
    print(f"{'host':<22}{'STRING':>9}{'TISSUES':>10}{'kept':>8}"
          f"{'DeepLoc':>10}{'kept':>8}{'both':>9}{'kept':>8}")
    for taxid, meta in hosts.items():
        total = len(pool[taxid])
        t, d = tis[taxid], dl[taxid]
        both = t & d
        print(f"{meta['label']:<22}{total:>9}{len(t):>10}{pct(len(t), total):>8}"
              f"{len(d):>10}{pct(len(d), total):>8}{len(both):>9}{pct(len(both), total):>8}")

    print("\n=== Human: TISSUES vs HPA single-cell expression ===")
    human = 9606
    total = len(pool[human])
    hpa_sets = hpa_pass(CONFIG, HPA_NTPM_CUTOFFS)
    hpa_any = hpa_any_tissue_universe(CONFIG)
    print(f"human STRING proteins                          {total:>7}")
    print(f"mappable to an HPA gene (any tissue)           {len(hpa_any):>7}  {pct(len(hpa_any), total)}")
    t = tis[human]
    print(f"pass TISSUES (>= {main.TISSUE_CUTOFF})                          {len(t):>7}  {pct(len(t), total)}")
    for c in HPA_NTPM_CUTOFFS:
        h = hpa_sets[c]
        print(f"pass HPA in a config tissue (nTPM > {c})       {len(h):>7}  {pct(len(h), total)}")

    print("\n=== Agreement on human (restricted to HPA-mappable proteins) ===")
    for c in HPA_NTPM_CUTOFFS:
        h = hpa_sets[c]
        tm, hm = t & hpa_any, h & hpa_any
        print(f"\n nTPM > {c}:")
        print(f"  both filters keep            {len(tm & hm):>7}")
        print(f"  TISSUES only (HPA drops)     {len(tm - hm):>7}")
        print(f"  HPA only (TISSUES drops)     {len(hm - tm):>7}")
        print(f"  neither                      {len(hpa_any - tm - hm):>7}")
        print(f"  of TISSUES-passing, HPA agrees {pct(len(tm & hm), len(tm))}")
        print(f"  of HPA-passing, TISSUES agrees {pct(len(tm & hm), len(hm))}")

    print("\n=== What the pipeline actually loses to the TISSUES filter, after DeepLoc ===")
    for taxid, meta in hosts.items():
        d = dl[taxid]
        both = d & tis[taxid]
        print(f"{meta['label']:<22} DeepLoc-surface {len(d):>6} -> {len(both):>6} "
              f"after TISSUES  ({pct(len(both), len(d))} kept)")

    print("\n=== Final pool at each TISSUES cutoff (after DeepLoc) ===")
    print(f"{'host':<22}" + ''.join(f"{'>=' + str(c):>10}" for c in SCORE_SWEEP))
    for taxid, meta in hosts.items():
        d = dl[taxid]
        counts = [len({p for p, v in scores[taxid].items() if v >= c} & d) for c in SCORE_SWEEP]
        print(f"{meta['label']:<22}" + ''.join(f"{n:>10}" for n in counts))

    print("\n=== HPA coverage of non-human hosts ===")
    print("HPA single-cell tissue data is human only; the other hosts have no HPA")
    print("equivalent in the config, so an HPA-only filter would leave them unfiltered.")


if __name__ == '__main__':
    main_(sys.argv[1] if len(sys.argv) > 1 else 'experiments')
