'''
AlphaFold 3 Server jobs for a benchmark of the confidence score against ipTM. Predicted
links are sampled evenly across confidence bins and, within each, between links whose
evidence is mostly experimental and those mostly from curated databases, one link per
pair of orthologous groups. Links spread over the strata are also run with a random
host protein of similar length, reachable in the same compartments and not predicted
to bind the parasite protein, as controls. Writes the jobs not yet run, predicted links
first, as JSON files of at most --batch-size jobs to upload one a day at
alphafoldserver.com, and a table of them with the batch of each. The batches are
fixed once written: a rerun reports how far each has run, and --rebatch splits the
jobs still to run anew.
'''
import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts import alphafold_results
from scripts.build_alphafold_jobs import controls, job, job_name, sequences

BINS = [0.35, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


def sample(predictions, seqs, per_stratum, max_tokens, rng):
    '''
    Links per confidence bin and evidence source, at most one per pair of orthologous
    groups and none longer than the server should be given.
    '''
    links = predictions[predictions['edge_type'] == 'inter-species'].copy()
    for column in ('weight', 'experimental_evidence_score', 'databases_evidence_score'):
        links[column] = links[column].astype(float)
    links['tokens'] = [len(seqs[s]) + len(seqs[t]) if s in seqs and t in seqs else np.inf
                       for s, t in zip(links['source'], links['target'])]
    links = links[links['tokens'] <= max_tokens]
    links['bin'] = pd.cut(links['weight'], BINS, right=False, include_lowest=True)
    links['evidence'] = np.where(links['experimental_evidence_score'] >=
                                 links['databases_evidence_score'], 'experimental', 'databases')
    # one link per group pair, so that orthologs of one link don't fill a stratum
    links = links.sample(frac=1, random_state=rng.integers(2**32)).drop_duplicates(['group1', 'group2'])
    picked = []
    for (_, _), stratum in links.groupby(['bin', 'evidence'], observed=True):
        # a parasite's nth link is taken only after every parasite's (n-1)th
        turn = stratum.groupby('taxid1').cumcount()
        picked.append(stratum.iloc[np.argsort(turn.to_numpy(), kind='stable')].head(per_stratum))
        if len(stratum) < per_stratum:
            print(f'{stratum["bin"].iloc[0]} {stratum["evidence"].iloc[0]}: '
                  f'only {len(stratum)} links')

    return pd.concat(picked)


def write_batches(table, jobs, order, out, batch_size, results, rebatch):
    '''
    Writes the jobs not yet run as JSON files of at most batch_size, in the order given by
    order(rows), and the table with the batch of each job to out.csv. The batches are fixed
    once written: a rerun keeps them, unless rebatch, and reports how far each has run.
    Batch 0 holds the jobs run before the batches were fixed, or run in another set.
    '''
    os.makedirs(os.path.dirname(out), exist_ok=True)
    table = table.copy()
    table['run'] = [alphafold_results.job_folder(n, results) is not None for n in table['job']]
    previous = pd.read_csv(f'{out}.csv') if os.path.exists(f'{out}.csv') else None
    if (not rebatch and previous is not None and 'batch' in previous.columns
            and list(previous['job']) == list(table['job'])):
        table['batch'] = previous['batch']
    else:
        pending = order(table[~table['run'] & table['job'].isin(jobs)])
        table['batch'] = 0
        table.loc[pending.index, 'batch'] = np.arange(len(pending)) // batch_size + 1
        for old in glob.glob(f'{out}_[0-9][0-9].json'):
            os.remove(old)
        for b, names in pending.groupby(table.loc[pending.index, 'batch'], sort=False)['job']:
            with open(f'{out}_{b:02d}.json', 'w') as f:
                json.dump([jobs[n] for n in names], f, indent=1)
    table.drop(columns='run').to_csv(f'{out}.csv', index=False)
    for b, batch in table[table['batch'] > 0].groupby('batch'):
        state = ('done' if batch['run'].all() else
                 f"{batch['run'].sum()} of {len(batch)} run" if batch['run'].any() else 'to run')
        kinds = ' and '.join(sorted(set(batch['kind'])))
        print(f'{out}_{b:02d}.json: {len(batch)} jobs, {kinds}: {state}')
    print(f"{table['run'].sum()} of {len(table)} jobs run")

    return table


def build(data_dir, per_stratum, n_controls, max_tokens, length_tolerance, batch_size,
          seed, results, rebatch, out):
    predictions = pd.read_parquet(os.path.join(data_dir, 'predictions.parquet'))
    eligible = pd.read_parquet(os.path.join(data_dir, 'eligible_proteins.parquet'))
    eligible['taxid'] = eligible['taxid'].astype(str)
    seqs = sequences(data_dir, set(predictions['taxid1']) | set(predictions['taxid2']))
    rng = np.random.default_rng(seed)
    links = sample(predictions, seqs, per_stratum, max_tokens, rng)

    # the controls go to links spread evenly over the strata
    with_control = set(np.linspace(0, len(links) - 1, min(n_controls, len(links))).round().astype(int))
    predicted, control_jobs, rows = [], [], []
    for i, link in enumerate(links.itertuples()):
        host = eligible[eligible['taxid'] == link.taxid2].set_index('protein')
        common = {'link': i, 'parasite': link.taxid1_label, 'host': link.taxid2_label,
                  'bin': str(link.bin), 'evidence': link.evidence, 'score': link.weight,
                  'experimental score': link.experimental_evidence_score,
                  'databases score': link.databases_evidence_score,
                  'group A': link.group1, 'group B': link.group2}
        partners = [('predicted', link.target, link.target_name, predicted)]
        if i in with_control:
            control = controls(link.target, link.source, host, predictions, seqs, 1,
                               length_tolerance, rng)[0]
            partners.append(('control', control, host.loc[control, 'name'], control_jobs))
        for kind, partner, partner_name, jobs in partners:
            name = job_name(f'bm{i:03d}', kind, link.source_name, partner_name)
            jobs.append(job(name, [link.source, partner], seqs))
            rows.append({'job': name, 'kind': kind, **common, 'chain A': link.source,
                         'chain B': partner, 'name A': link.source_name,
                         'name B': partner_name,
                         'tokens': len(seqs[link.source]) + len(seqs[partner])})

    # predicted links first, so their results come in before the controls', and each
    # stratum's nth job after every stratum's (n-1)th, so that every batch spans them all
    def order(pending):
        pending = pending.assign(first=pending['kind'] != 'predicted',
                                 turn=pending.groupby(['kind', 'bin', 'evidence']).cumcount())
        return pending.sort_values(['first', 'turn', 'link'])

    table = write_batches(pd.DataFrame(rows), {j['name']: j for j in predicted + control_jobs},
                          order, out, batch_size, results, rebatch)
    print(f'{len(links)} links and {len(control_jobs)} controls')
    print(table[table['kind'] == 'predicted'].groupby(['bin', 'evidence']).size().unstack())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--per-stratum', type=int, default=8,
                        help='links per confidence bin and evidence source')
    parser.add_argument('--controls', type=int, default=30,
                        help='links run against a random host protein too')
    parser.add_argument('--max-tokens', type=int, default=1500,
                        help='longest pair, in residues, to sample')
    parser.add_argument('--length-tolerance', type=float, default=0.1,
                        help='fraction a control may differ in length from the partner')
    parser.add_argument('--batch-size', type=int, default=30,
                        help='jobs per JSON file, the number the server runs in a day')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--results', default=alphafold_results.RESULTS,
                        help='where the server downloads are unpacked; jobs found there are '
                             'left out of new batches')
    parser.add_argument('--rebatch', action='store_true',
                        help='split the jobs not yet run into new batch files; otherwise '
                             'the batches written before are kept')
    parser.add_argument('--out', default='paper/alphafold/benchmark',
                        help='output path without extension')
    args = parser.parse_args()
    build(args.data_dir, args.per_stratum, args.controls, args.max_tokens, args.length_tolerance,
          args.batch_size, args.seed, args.results, args.rebatch, args.out)
