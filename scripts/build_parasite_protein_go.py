'''
The parasite processes behind the predicted interactions. For each parasite, the GO
biological processes over-represented among its proteins with a predicted host
interaction are tested against the parasite's own proteins the pipeline passed
(web_utils.go_background), with the one-sided Fisher's exact test and Benjamini-Hochberg
correction the app uses. Each process also counts the parasite orthologous groups its
proteins come from, as one expanded family can carry a process on its own.
Writes paper/tables/parasite_protein_go_enriched.csv (every process enriched for every
parasite).
'''
import argparse
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [ROOT, os.path.join(ROOT, 'app')]

import utils
import web_utils


def enrich(config, data_dir, edges):
    '''Every process tested for every parasite of these links, one row per parasite and term.'''
    tests = []
    for taxid, rows in edges.groupby('taxid1'):
        test = utils.calculate_enrichment(set(rows['source']),
                                          web_utils.go_background(data_dir, (taxid,)))
        if not test.empty:
            tests.append(test.drop(columns=['B', 'C', 'D'])
                         .rename(columns={'A': 'n_proteins'})
                         .assign(taxid=int(taxid), niche=web_utils.parasite_niche(config, taxid)))

    return pd.concat(tests, ignore_index=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--host', type=int, default=9606, help='host taxid (default: human)')
    parser.add_argument('--min-score', type=float, default=0.35,
                        help='links below this confidence score are left out')
    parser.add_argument('--fdr', type=float, default=0.05)
    parser.add_argument('--table', default='paper/tables/parasite_protein_go_enriched.csv')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    parasites = config['parasites']
    predictions = pd.read_parquet(os.path.join(args.data_dir, 'predictions.parquet'))
    edges = predictions[(predictions['taxid2'] == str(args.host))
                        & (predictions['weight'].astype(float) >= args.min_score)]
    enriched = enrich(config, args.data_dir, edges)
    enriched = enriched[enriched['fdr_bh'] <= args.fdr].copy()

    group_of = edges.drop_duplicates('source').set_index('source')['group1']
    enriched['n_groups'] = enriched['nodes'].map(lambda nodes: group_of[nodes.split(',')].nunique())
    enriched['parasite'] = enriched['taxid'].map(lambda t: parasites[t]['label'])
    enriched = enriched.sort_values(['niche', 'parasite', 'fdr_bh'])
    enriched[['parasite', 'taxid', 'niche', 'go_id', 'go_term', 'n_proteins', 'n_groups',
              'odds_ratio', 'p_value', 'fdr_bh', 'nodes']].to_csv(args.table, index=False)
    print(f"Wrote {args.table}: {len(enriched)} enriched processes over "
          f"{enriched['taxid'].nunique()} parasites, {enriched['go_id'].nunique()} distinct")
