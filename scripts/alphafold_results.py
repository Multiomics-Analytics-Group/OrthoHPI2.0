'''
Reads AlphaFold 3 Server results unpacked under data/alphafold/folds_*/: each job's
models with their summary confidences and the ipSAE of their interface.
'''
import glob
import json
import os

import numpy as np

RESULTS = os.path.join('data', 'alphafold')
# PAE, in Å, below which a residue pair counts toward ipSAE (Dunbrack 2025)
IPSAE_CUTOFF = 10


def job_folder(job, results=RESULTS):
    '''The folder of a job, in whichever download it came in; None if not yet run.'''
    found = glob.glob(os.path.join(results, 'folds_*', job.lower().replace('-', '_')))

    return sorted(found)[-1] if found else None


def ipsae(pae, chains, cutoff=IPSAE_CUTOFF):
    '''
    ipSAE of a two-chain model (Dunbrack 2025): for each residue of one chain, a TM-score
    over the residues of the other chain it is aligned to within the cutoff, with d0 from
    their number; the best residue, in the better of the two directions.
    '''
    best = 0.0
    for a in np.unique(chains):
        rows = pae[np.ix_(chains == a, chains != a)]
        for row in rows:
            close = row[row < cutoff]
            if close.size:
                d0 = max(1.0, 1.24 * (max(close.size, 27) - 15) ** (1 / 3) - 1.8)
                best = max(best, float(np.mean(1 / (1 + (close / d0) ** 2))))

    return best


def models(folder):
    '''The confidences of a job's models, best-ranked first.'''
    rows = []
    for rank, path in enumerate(sorted(glob.glob(os.path.join(folder, '*_summary_confidences_*.json')))):
        summary = json.load(open(path))
        full = json.load(open(path.replace('summary_confidences', 'full_data')))
        rows.append({'model': rank, 'iptm': summary['iptm'], 'ptm': summary['ptm'],
                     'ipsae': ipsae(np.array(full['pae']), np.array(full['token_chain_ids'])),
                     'ranking score': summary['ranking_score'],
                     'fraction disordered': summary['fraction_disordered'],
                     'min PAE A-B': summary['chain_pair_pae_min'][0][1],
                     'min PAE B-A': summary['chain_pair_pae_min'][1][0]})

    return rows
