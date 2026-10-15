'''
How strong the evidence of direct binding is for each human reference of the
direct-interaction benchmark (scripts/build_alphafold_direct_jobs.py), and whether the
comparison of predicted and reference models holds when only strong evidence counts.
Strong is IntAct's direct interaction or BioGRID's co-crystal structure, far western or
reconstituted complex; suggestive is two-hybrid, protein-fragment complementation,
protein-peptide, biochemical activity or an IntAct enzymatic reaction, cleavage or
covalent bond; FRET alone is proximity, not binding. Writes
paper/tables/alphafold_direct_evidence.csv.
'''
import argparse
import os
import sys

import pandas as pd
from scipy.stats import wilcoxon

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.binding_evidence import interactions

STRONG = {'direct interaction', 'Co-crystal Structure', 'Far Western', 'Reconstituted Complex'}
PROXIMITY = {'FRET'}


def strength(evidence):
    if evidence & STRONG:
        return 'strong'
    if evidence - PROXIMITY:
        return 'suggestive'

    return 'proximity only' if evidence else 'none'


def compare(links, label):
    for metric in ('iptm', 'ipsae'):
        reference, predicted = links[f'{metric} reference'], links[f'{metric} predicted']
        print(f'{label}, {metric}: n = {len(links)}, reference higher in '
              f'{(reference > predicted).sum()}, predicted in {(predicted > reference).sum()}; '
              f'median {reference.median():.2f} vs {predicted.median():.2f}, '
              f'Wilcoxon p = {wilcoxon(reference, predicted).pvalue:.2g}')


def build(jobs_file, scores_file, data_dir, out):
    known, releases = interactions(os.path.join(data_dir, 'downloads'))
    direct = known[known['kind'] == 'direct']
    # IntAct's evidence by its interaction type, BioGRID's by its experimental system
    direct = direct.assign(evidence=direct['type'].where(direct['source'] == 'IntAct',
                                                         direct['method']))
    by_pair = direct.groupby('pair')['evidence'].agg(set)
    references = pd.read_csv(jobs_file)
    references = references[references['kind'] == 'reference']
    rows = []
    for ref in references.to_dict('records'):
        pair = '|'.join(sorted((ref['name A'].upper(), ref['name B'].upper())))
        evidence = by_pair.get(pair, set())
        rows.append({'link': ref['link'], 'name A': ref['name A'], 'name B': ref['name B'],
                     'evidence': '; '.join(sorted(evidence)), 'strength': strength(evidence)})
    table = pd.DataFrame(rows).merge(pd.read_csv(scores_file)[
        ['link', 'same partner', 'iptm predicted', 'iptm reference', 'ipsae predicted',
         'ipsae reference']], on='link')
    table.assign(releases=releases).to_csv(out, index=False)

    print(table['strength'].value_counts().to_string())
    strong = table[table['strength'] == 'strong']
    compare(strong, 'strong')
    compare(strong[strong['same partner']], 'strong, same partner')
    compare(table[table['strength'] != 'proximity only'], 'without FRET-only')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--jobs', default='paper/alphafold/direct.csv')
    parser.add_argument('--scores', default='paper/tables/alphafold_direct.csv',
                        help='the interface confidences of scripts/build_alphafold_direct_figure.py')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--out', default='paper/tables/alphafold_direct_evidence.csv')
    args = parser.parse_args()
    build(args.jobs, args.scores, args.data_dir, args.out)
