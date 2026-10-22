'''
The options for a case study of one parasite in one host, narrowed as a funnel: every
predicted link, the links of one tissue the parasite infects, those whose host protein is
enhanced in one cell type of that tissue, and those whose host protein is found in none
of the parasite's other tissues. For each tissue and cell type it counts the links and
proteins at every stage and lists the host processes most enriched there, so a tissue and
cell type can be picked for the figure (scripts/build_pf_case_study.py). Writes paper/tables/pf_case_study_options.csv.

It also reports where the experimentally characterized interactions of the pair stand in
the pipeline (KNOWN_PAIRS): whether each side passed the filters and joined an orthology
group, and whether the link was predicted.
'''
import argparse
import gzip
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [ROOT, os.path.join(ROOT, 'app')]

import utils
import web_utils

# experimentally characterized P. falciparum-human interactions, as (parasite UniProt,
# host UniProt, label, source); the BioGRID MV release is 5.0.261, IntAct as of 2026-09
KNOWN_PAIRS = [
    ('Q8IFM5', 'P35613', 'RH5-BSG', 'IntAct, PMID 24297912'),
    ('Q8IBE8', 'P02724', 'EBA-175-GYPA', 'IntAct, PMID 24297912'),
    ('Q8II38', 'Q9BZL1', 'Prp38-UBL5', 'BioGRID MV, PMID 21614000'),
]
TOP_TERMS = 5
# HPA's specificity rule: a gene is enhanced in a cell type where it is detected and its
# nCPM is at least this many times the mean over the other cell types, here of one tissue
ENHANCED_FOLD = 4
DETECTED = 1.0


def enhanced(annotation, tissue):
    '''The host proteins enhanced in each cell type of one tissue, as {cell type: set}.'''
    cells = annotation[annotation['Tissue'] == tissue].dropna(subset=['Cell type'])
    ncpm = cells.groupby(['Gene', 'Cell type'])['nCPM'].max().unstack().fillna(0)
    others = (ncpm.sum(axis=1).to_numpy()[:, None] - ncpm) / (ncpm.shape[1] - 1)
    calls = (ncpm > DETECTED) & (ncpm >= ENHANCED_FOLD * others)

    return {cell_type: set(calls.index[calls[cell_type]]) for cell_type in calls.columns}


def top_terms(proteins, go_df):
    '''The host processes most enriched among these proteins, by FDR.'''
    enrichment = utils.calculate_enrichment(set(proteins), go_df)
    if enrichment.empty:
        return ''
    significant = enrichment[enrichment['fdr_bh'] < 0.05].head(TOP_TERMS)

    return '; '.join(f"{row.go_term} ({row.fdr_bh:.1e})" for row in significant.itertuples())


def stage(label, tissue, cell_type, edges, go_df):
    return {'stage': label, 'tissue': tissue, 'cell_type': cell_type, 'links': len(edges),
            'parasite_proteins': edges['source'].nunique(),
            'host_proteins': edges['target'].nunique(),
            'top_go': top_terms(edges['target'], go_df)}


def load(config, data_dir, parasite, host):
    '''
    The predicted links of the pair, the tissue annotation of their host proteins, and
    the host GO annotation of the pool the parasite's niche reaches, which is the
    background the app tests enrichment against.
    '''
    predictions = pd.read_parquet(os.path.join(data_dir, 'predictions.parquet'))
    edges = predictions[(predictions['taxid1'] == str(parasite))
                        & (predictions['taxid2'] == str(host))].drop_duplicates(['source', 'target'])
    annotation = pd.read_parquet(os.path.join(data_dir, 'tissues_cell_types.parquet'))
    annotation = annotation[annotation['Gene'].isin(edges['target'])]
    go_df = pd.read_parquet(os.path.join(data_dir, 'gos.parquet'), filters=[('taxid', '==', host)])
    pool = web_utils.filtered_pool(data_dir, (str(host),),
                                   niche=web_utils.parasite_niche(config, parasite))

    return edges, annotation, go_df[go_df['#string_protein_id'].isin(pool)]


def tissue_proteins(config, annotation, parasite):
    '''The host proteins of each tissue the parasite infects, in the order config gives.'''
    tissues = [config['tissues'][t] for t in config['parasites'][parasite]['tissues']]

    return {t: set(annotation.loc[annotation['Tissue'] == t, 'Gene']) for t in tissues}


def build(config, data_dir, parasite, host, min_host_proteins):
    edges, annotation, go_df = load(config, data_dir, parasite, host)
    in_tissue = tissue_proteins(config, annotation, parasite)
    tissues = list(in_tissue)

    rows = [stage('network', '', '', edges, go_df)]
    for tissue in tissues:
        others = set().union(*(in_tissue[t] for t in tissues if t != tissue))
        tissue_edges = edges[edges['target'].isin(in_tissue[tissue])]
        specific = tissue_edges[~tissue_edges['target'].isin(others)]
        rows.append(stage('tissue', tissue, '', tissue_edges, go_df))
        rows.append(stage('tissue-specific', tissue, '', specific, go_df))
        for cell_type, targets in enhanced(annotation, tissue).items():
            if len(targets) < min_host_proteins:
                continue
            cell_edges = tissue_edges[tissue_edges['target'].isin(targets)]
            rows.append(stage('cell type', tissue, cell_type, cell_edges, go_df))
            rows.append(stage('cell type, tissue-specific', tissue, cell_type,
                              cell_edges[~cell_edges['target'].isin(others)], go_df))

    return pd.DataFrame(rows), edges


def known_pairs(data_dir, parasite, edges):
    '''Where each characterized pair stands: filters passed, orthologous group, prediction.'''
    eligible = pd.read_parquet(os.path.join(data_dir, 'eligible_proteins.parquet'))
    uniprot = pd.concat([edges[['source', 'source_uniprot']].set_axis(['id', 'uniprot'], axis=1),
                         edges[['target', 'target_uniprot']].set_axis(['id', 'uniprot'], axis=1)])
    ids = dict(zip(uniprot['uniprot'], uniprot['id']))
    members = {}
    with gzip.open(os.path.join(data_dir, 'downloads', '2759_members.tsv.gz'), 'rt') as handle:
        wanted = {f'{parasite}.{p}' for p, _, _, _ in KNOWN_PAIRS}
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            for protein in wanted.intersection(fields[4].split(',')):
                members[protein] = fields[1]
    rows = []
    for parasite_acc, host_acc, label, source in KNOWN_PAIRS:
        parasite_id = f'{parasite}.{parasite_acc}'
        predicted = ((edges['source_uniprot'] == parasite_acc)
                     & (edges['target_uniprot'] == host_acc)).any()
        rows.append({'pair': label, 'source': source,
                     'parasite_eligible': parasite_id in set(eligible['protein']),
                     'parasite_group': members.get(parasite_id, ''),
                     'host_in_network': host_acc in ids, 'predicted': predicted})

    return pd.DataFrame(rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--parasite', type=int, default=36329,
                        help='parasite taxid (default: Plasmodium falciparum)')
    parser.add_argument('--host', type=int, default=9606, help='host taxid (default: human)')
    parser.add_argument('--min-host-proteins', type=int, default=5,
                        help='cell types enhancing fewer of the host proteins are left out')
    parser.add_argument('--output', default='paper/tables/pf_case_study_options.csv')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    options, edges = build(config, args.data_dir, args.parasite, args.host, args.min_host_proteins)
    options.to_csv(args.output, index=False)
    with pd.option_context('display.width', 250, 'display.max_colwidth', 120,
                           'display.max_rows', 500):
        print(options.to_string(index=False))
        print()
        print(known_pairs(args.data_dir, args.parasite, edges).to_string(index=False))
    print(f'Wrote {args.output}')
