'''
AlphaFold 3 Server jobs for a benchmark of predicted human links whose
orthologous groups have human members shown to bind directly in BioGRID or IntAct
(scripts/binding_evidence.py). Links are sampled evenly across confidence bins,
one per pair of orthologous groups, from the parasites in turn. Each is run with its
reference: the human pair of its groups with direct evidence, preferring the predicted
partner itself, then the pair with the most publications. Direct human links of an
earlier benchmark of links sampled regardless of evidence, listed in
paper/alphafold/benchmark.csv, are taken over with their models, so they were not run
again. Writes the jobs not yet run as JSON files of at most
--batch-size jobs, a link and its reference together, and a table of them with the
batch of each; the batches are fixed once written.
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
from scripts.binding_evidence import HUMAN, human_members, interactions
from scripts.build_alphafold_jobs import job, job_name, sequences

# the confidence-score bins links are sampled across
BINS = [0.35, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


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


def direct_pairs(links, members, known, results):
    '''
    For each group pair, its human pairs with direct evidence and their publications,
    the predicted partner's first, then the best-published, then one already run.
    '''
    direct = known[known['kind'] == 'direct']
    publications = direct.groupby('pair')['publication'].nunique()
    found = {}
    for a, b, target in zip(links['group1'], links['group2'], links['target_name']):
        if (a, b) in found:
            continue
        pairs = [(x, y, publications['|'.join(sorted((x, y)))])
                 for x in members.get(a, ()) for y in members.get(b, ())
                 if x != y and '|'.join(sorted((x, y))) in publications.index]
        # among equals, a reference already run comes first, so a rerun keeps the models
        # there are; then by name, as the members come in sets, in no fixed order
        found[(a, b)] = sorted(pairs, key=lambda p: (
            p[1] != target.upper(), -p[2],
            alphafold_results.job_folder(f'dr*_reference_{p[0]}_{p[1]}', results) is None,
            p[0], p[1]))

    return found


def build(data_dir, per_bin, max_tokens, max_reference_tokens, batch_size, seed, results,
          rebatch, benchmark, out):
    downloads = os.path.join(data_dir, 'downloads')
    predictions = pd.read_parquet(os.path.join(data_dir, 'predictions.parquet'))
    links = predictions[(predictions['edge_type'] == 'inter-species') &
                        (predictions['taxid2'] == str(HUMAN))].copy()
    links['weight'] = links['weight'].astype(float)
    eligible = pd.read_parquet(os.path.join(data_dir, 'eligible_proteins.parquet'))
    human = eligible[eligible['taxid'].astype(str) == str(HUMAN)]
    symbols = dict(zip(human['protein'], human['name']))
    by_symbol = {name.upper(): protein for protein, name in symbols.items()}
    members = human_members(downloads, set(links['group1']) | set(links['group2']), symbols)
    known, releases = interactions(downloads)
    pairs = direct_pairs(links, members, known, results)
    seqs = sequences(data_dir, set(links['taxid1']) | {str(HUMAN)})

    # a link needs a reference; human proteins shown to bind directly are often far longer
    # than the parasite protein, so references have their own limit
    def reference(link):
        for a, b, n in pairs[(link.group1, link.group2)]:
            chains = [by_symbol.get(a), by_symbol.get(b)]
            if (all(c in seqs for c in chains) and
                    sum(len(seqs[c]) for c in chains) <= max_reference_tokens):
                return a, b, n, chains

        return None

    links = links[[bool(pairs.get((a, b))) for a, b in zip(links['group1'], links['group2'])]]
    links = links[[s in seqs and t in seqs and len(seqs[s]) + len(seqs[t]) <= max_tokens
                   for s, t in zip(links['source'], links['target'])]]
    links['bin'] = pd.cut(links['weight'], BINS, right=False, include_lowest=True)

    # the direct human links of the first benchmark come first, with their jobs
    first = pd.read_csv(benchmark)
    first = first[(first['kind'] == 'predicted') & (first['host'] == 'Homo sapiens (human)')]
    taken = dict(zip(zip(first['chain A'], first['chain B']), first['job']))
    links['reused'] = [taken.get((s, t)) for s, t in zip(links['source'], links['target'])]

    rng = np.random.default_rng(seed)
    links = links.sample(frac=1, random_state=rng.integers(2**32))
    links = links.sort_values('reused', na_position='last', kind='stable')
    links = links.drop_duplicates(['group1', 'group2'])
    picked = []
    for _, in_bin in links.groupby('bin', observed=True):
        in_bin = in_bin.assign(ref=[reference(link) for link in in_bin.itertuples()])
        in_bin = in_bin[in_bin['ref'].notna()]
        # reused links first, then a parasite's nth link after every parasite's (n-1)th
        turn = np.where(in_bin['reused'].notna(), -1, in_bin.groupby('taxid1').cumcount())
        picked.append(in_bin.iloc[np.argsort(turn, kind='stable')].head(per_bin))
        if len(in_bin) < per_bin:
            print(f'{in_bin["bin"].iloc[0]}: only {len(in_bin)} links')
    links = pd.concat(picked)

    jobs, rows = {}, []
    for i, link in enumerate(links.itertuples()):
        a, b, n, ref_chains = link.ref
        common = {'link': i, 'parasite': link.taxid1_label, 'bin': str(link.bin),
                  'score': link.weight, 'group A': link.group1, 'group B': link.group2,
                  'reference publications': n, 'same partner': b == link.target_name.upper()}
        name = link.reused or job_name(f'dr{i:03d}', 'predicted', link.source_name,
                                       link.target_name)
        for kind, chains, labels in (
                ('predicted', [link.source, link.target], [link.source_name, link.target_name]),
                ('reference', ref_chains, [symbols[c] for c in ref_chains])):
            if kind == 'reference':
                name = job_name(f'dr{i:03d}', kind, *labels)
            if not (kind == 'predicted' and link.reused):
                jobs[name] = job(name, chains, seqs)
            rows.append({'job': name, 'kind': kind, **common, 'chain A': chains[0],
                         'chain B': chains[1], 'name A': labels[0], 'name B': labels[1],
                         'tokens': sum(len(seqs[c]) for c in chains),
                         'reused': kind == 'predicted' and bool(link.reused)})

    # each bin's nth link after every bin's (n-1)th, with its reference beside it
    def order(pending):
        turn = pending.groupby('bin')['link'].rank(method='dense')
        return pending.assign(turn=turn).sort_values(['turn', 'bin', 'link', 'kind'])

    table = write_batches(pd.DataFrame(rows).assign(releases=releases), jobs, order, out,
                          batch_size, results, rebatch)
    predicted = table[table['kind'] == 'predicted']
    print(f"{len(links)} links ({predicted['reused'].sum()} from the first benchmark) "
          f"from {links['taxid1'].nunique()} parasites; "
          f"{predicted['same partner'].sum()} references share the predicted partner")
    print(predicted.groupby('bin').size().to_string())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--per-bin', type=int, default=8, help='links per confidence bin')
    parser.add_argument('--max-tokens', type=int, default=1500,
                        help='longest predicted pair, in residues')
    parser.add_argument('--max-reference-tokens', type=int, default=2500,
                        help='longest reference pair, in residues')
    parser.add_argument('--batch-size', type=int, default=30,
                        help='jobs per JSON file, the number the server runs in a day')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--results', default=alphafold_results.RESULTS,
                        help='where the server downloads are unpacked; jobs found there are '
                             'left out of new batches')
    parser.add_argument('--rebatch', action='store_true',
                        help='split the jobs not yet run into new batch files; otherwise '
                             'the batches written before are kept')
    parser.add_argument('--benchmark', default='paper/alphafold/benchmark.csv',
                        help='the earlier benchmark, whose direct human links are taken over')
    parser.add_argument('--out', default='paper/alphafold/direct',
                        help='output path without extension')
    args = parser.parse_args()
    build(args.data_dir, args.per_bin, args.max_tokens, args.max_reference_tokens,
          args.batch_size, args.seed,
          args.results, args.rebatch, args.benchmark, args.out)
