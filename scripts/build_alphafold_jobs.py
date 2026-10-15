'''
AlphaFold 3 Server jobs for case-study links. Each predicted parasite-host pair is run
alongside the host-host complex it was transferred from, the parasite's host ortholog
taking its place, as a reference for the confidence and interface of a
real complex; and the parasite protein against random host proteins of similar length
reachable in the same compartments and not predicted to bind it, as negative controls.
Writes the jobs as one JSON to upload at alphafoldserver.com and a table of them.
'''
import argparse
import json
import os
import re

from Bio import SeqIO
import numpy as np
import pandas as pd

# (case study, parasite protein, predicted host partner, host ortholog of the parasite
# protein); the host-host references have solved structures: LC3B with ATG4B (PDB 2Z0D),
# the NBR1 LIR (2L8J) and the ATG3 LIR, PPP2R3B (PR70) with the PP2A A subunit (4I5L)
PAIRS = [
    ('lb', '5660.A4H5J3', 'ATG4B', 'MAP1LC3B'),
    ('lb', '5660.A4H5J3', 'NBR1', 'MAP1LC3B'),
    ('lb', '5660.A4H5J3', 'ATG3', 'MAP1LC3B'),
    ('pf', '36329.Q8IDD8', 'PPP2R1B', 'PPP2R3B'),
    ('pf', '36329.Q8I5V5', 'CIAO2A', 'CIAO1'),
]
HOST = '9606'
# tokens the server takes in one job
MAX_TOKENS = 5000


def sequences(data_dir, taxids):
    seqs = {}
    for taxid in taxids:
        path = os.path.join(data_dir, 'deeploc', 'input', f'{taxid}.fasta')
        seqs.update({r.id: str(r.seq) for r in SeqIO.parse(path, 'fasta')})
    return seqs


def controls(partner, parasite, eligible, predictions, seqs, n, length_tolerance, rng):
    '''
    Host proteins of the partner's length, within the tolerance, reachable in the same
    compartments, and not predicted to bind the parasite protein.
    '''
    row = eligible.loc[partner]
    length = len(seqs[partner])
    bound = set(predictions.loc[predictions['source'] == parasite, 'target'])
    pool = eligible[(eligible['reachable_extracellular'] == row['reachable_extracellular']) &
                    (eligible['reachable_intracellular'] == row['reachable_intracellular'])]
    pool = [p for p in pool.index if p in seqs and p not in bound and p != partner and
            abs(len(seqs[p]) - length) <= length_tolerance * length]
    return list(rng.choice(sorted(pool), size=n, replace=False))


def job_name(*parts):
    '''The parts joined into a name the server takes: letters, digits, _ and -.'''
    return re.sub(r'[^A-Za-z0-9_-]', '-', '_'.join(parts))


def job(name, chains, seqs):
    return {'name': name, 'modelSeeds': [], 'dialect': 'alphafoldserver', 'version': 1,
            'sequences': [{'proteinChain': {'sequence': seqs[c], 'count': 1}}
                          for c in chains]}


def build(data_dir, n_controls, length_tolerance, seed, out):
    predictions = pd.read_parquet(os.path.join(data_dir, 'predictions.parquet'))
    eligible = pd.read_parquet(os.path.join(data_dir, 'eligible_proteins.parquet'))
    eligible = eligible[eligible['taxid'].astype(str) == HOST].set_index('protein')
    by_name = pd.Series(eligible.index, index=eligible['name'])
    seqs = sequences(data_dir, {p.split('.')[0] for _, p, _, _ in PAIRS} | {HOST})
    names = dict(zip(predictions['source'], predictions['source_name']))
    rng = np.random.default_rng(seed)

    jobs, rows = [], []

    def add(study, pair, kind, chains, labels, score=None):
        name = job_name(study, kind, *labels)
        tokens = sum(len(seqs[c]) for c in chains)
        if tokens > MAX_TOKENS:
            raise SystemExit(f'{name} has {tokens} tokens, more than the server takes')
        jobs.append(job(name, chains, seqs))
        rows.append({'job': name, 'case study': study, 'pair': pair, 'kind': kind,
                     'chain A': chains[0], 'chain B': chains[1], 'name A': labels[0],
                     'name B': labels[1], 'tokens': tokens, 'score': score})

    for study, parasite, partner_name, ortholog_name in PAIRS:
        partner, ortholog = by_name[partner_name], by_name[ortholog_name]
        link = predictions[(predictions['source'] == parasite) &
                           (predictions['target'] == partner)]
        if link.empty:
            raise SystemExit(f'no predicted link between {parasite} and {partner_name}')
        label = names[parasite]
        pair = f'{label} – {partner_name}'
        add(study, pair, 'predicted', [parasite, partner], [label, partner_name],
            float(link['weight'].iloc[0]))
        add(study, pair, 'reference', [ortholog, partner], [ortholog_name, partner_name])
        for control in controls(partner, parasite, eligible, predictions, seqs,
                                n_controls, length_tolerance, rng):
            add(study, pair, 'control', [parasite, control], [label, eligible.loc[control, 'name']])

    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(f'{out}.json', 'w') as f:
        json.dump(jobs, f, indent=1)
    pd.DataFrame(rows).to_csv(f'{out}.csv', index=False)
    print(f'{len(jobs)} jobs written to {out}.json')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--controls', type=int, default=1,
                        help='random host proteins run against each parasite protein per pair')
    parser.add_argument('--length-tolerance', type=float, default=0.1,
                        help='fraction a control may differ in length from the partner')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--out', default='paper/alphafold/jobs',
                        help='output path without extension')
    args = parser.parse_args()
    build(args.data_dir, args.controls, args.length_tolerance, args.seed, args.out)
