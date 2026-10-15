'''
The depth of the paired alignments the AlphaFold 3 Server built for the
direct-interaction benchmark (scripts/build_alphafold_direct_jobs.py): for each job,
the rows of its paired alignment, and those in which both chains have a homolog rather
than a row of gaps, compared between each predicted link and its human reference and
against the interface confidence of their models
(scripts/build_alphafold_direct_figure.py). Writes
paper/tables/alphafold_direct_msa_depth.csv.
'''
import argparse
import glob
import os
import sys

import pandas as pd
from scipy.stats import spearmanr, wilcoxon

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts import alphafold_results

# fraction of a row's aligned positions that must be residues for it to count as a homolog
MIN_RESIDUES = 0.1


def a3m(path):
    '''The sequences of an A3M alignment, in order.'''
    sequences, current = [], []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if current:
                    sequences.append(''.join(current))
                current = []
            else:
                current.append(line.strip())
    if current:
        sequences.append(''.join(current))

    return sequences


def homolog(row):
    # in A3M, upper-case letters and gaps are the aligned positions, lower case insertions
    aligned = sum(c.isupper() or c == '-' for c in row)
    return sum(c.isupper() for c in row) >= MIN_RESIDUES * max(aligned, 1)


def depth(folder):
    '''The rows of a job's paired alignment, and those with a homolog of both chains.'''
    files = [glob.glob(os.path.join(folder, 'msas', f'*_paired_msa_chains_{c}.a3m'))
             for c in 'ab']
    if not all(files):
        return None, None
    a, b = a3m(files[0][0]), a3m(files[1][0])

    return len(a), sum(homolog(x) and homolog(y) for x, y in zip(a, b))


def build(jobs_file, scores_file, results, out):
    jobs = pd.read_csv(jobs_file)
    rows = []
    for job in jobs.to_dict('records'):
        folder = alphafold_results.job_folder(job['job'], results)
        if folder is not None:
            total, paired = depth(folder)
            rows.append({'link': job['link'], 'kind': job['kind'], 'rows': total,
                         'paired': paired})
    wide = pd.DataFrame(rows).pivot(index='link', columns='kind', values=['rows', 'paired'])
    wide.columns = [f'{measure} {kind}' for measure, kind in wide.columns]
    wide = wide.dropna().astype(int).reset_index()
    scores = pd.read_csv(scores_file)[['link', 'score', 'same partner', 'ipsae predicted',
                                       'ipsae reference']]
    table = wide.merge(scores, on='link')
    table.to_csv(out, index=False)

    deeper = (table['paired reference'] > table['paired predicted']).sum()
    print(f"{len(table)} links; median paired rows {table['paired predicted'].median():.0f} "
          f"predicted vs {table['paired reference'].median():.0f} reference; reference "
          f"deeper in {deeper}, Wilcoxon p = "
          f"{wilcoxon(table['paired reference'], table['paired predicted']).pvalue:.2g}")
    for kind in ('predicted', 'reference'):
        rho, p = spearmanr(table[f'paired {kind}'], table[f'ipsae {kind}'])
        print(f'{kind}: ipSAE against paired rows, rho = {rho:.2f} (p = {p:.2g})')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--jobs', default='paper/alphafold/direct.csv')
    parser.add_argument('--scores', default='paper/tables/alphafold_direct.csv',
                        help='the interface confidences of scripts/build_alphafold_direct_figure.py')
    parser.add_argument('--results', default=alphafold_results.RESULTS,
                        help='where the server downloads are unpacked, as folds_*')
    parser.add_argument('--out', default='paper/tables/alphafold_direct_msa_depth.csv')
    args = parser.parse_args()
    build(args.jobs, args.scores, args.results, args.out)
