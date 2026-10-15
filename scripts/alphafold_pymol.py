'''
PyMOL commands for the AlphaFold 3 Server models under data/alphafold/folds_*/. Run
from the repository root, `pymol scripts/alphafold_pymol.py`, or `run
scripts/alphafold_pymol.py` inside PyMOL, opens a chosen set of complexes, one scene each
(Page Up/Down to step through them), and adds the commands below. The set's objects are
numbered by complex, a predicted model sharing its number with its reference and,
where one is solved, the experimental structure of the host-host complex, e.g.
04_pred_atg8c_1_atg4b, 04_ref_map1lc3b_atg4b and 04_pdb_2z0d. Fetched entries are kept
in data/alphafold/pdb/.

    af_list [pattern]               jobs with models, e.g. af_list ciao
    af_load job [, model [, min_plddt]]
                                    a model, parasite protein blue, host protein grey,
                                    interface residues as sticks, residues of pLDDT
                                    below 50 (or min_plddt) hidden
    af_compare predicted, reference a parasite-host model and its host-host reference,
                                    the reference's host ortholog orange, superposed on
                                    the host partner they share (chain B) where the
                                    reference binds it
    af_pdb entry, chain, target [, others]
                                    a PDB entry's first model, its chain superposed on
                                    the target selection, grey, with its other chains
                                    given (e.g. B+C) violet, e.g.
                                    af_pdb 2z0d, A, 04_pred_atg8c_1_atg4b and chain B, B
    af_plddt [selection]            colour by pLDDT, as the AlphaFold database does

To fit an object to the window, `zoom <object> and visible`: plain `zoom` also fits the
hidden residues. The commands take the result folders' names, e.g. pf_predicted_pf3d7_1209400_ciao2a.
'''
import glob
import os

from pymol import CmdException, cmd, util

# PyMOL's __file__ is its own package, so the results are found from where it was started
RESULTS = os.path.join(os.getcwd(), 'data', 'alphafold')

PDB = os.path.join(RESULTS, 'pdb')

PARASITE = '0x2a78d6'
HOST = '0xb4b2a9'
ORTHOLOG = '0xeb6834'
EXPERIMENTAL = '0x4a3aa7'
# Å between heavy atoms of the two chains for a residue to count as interface
INTERFACE = 5.0
# pLDDT below which residues are hidden, and above which they are used to superpose
MIN_PLDDT = 50
CONFIDENT = 70
# Å from the reference's host ortholog within which the partner is superposed
SITE = 12

# the kinds of job, as object names shorten them
KINDS = {'predicted': 'pred', 'reference': 'ref', 'control': 'ctrl'}
# (predicted, reference or None, experimental structure or None) opened by `pymol
# scripts/alphafold_pymol.py`, numbered in this order; the benchmark links have no
# reference. A structure is (PDB entry, its chain superposed, the predicted model's chain
# it is superposed on, its other chains kept, the one touching the first if several):
# 4I5L PP2A Aα (A) with PR70 (B); 2L8J GABARAPL1 (A) with the NBR1 LIR (B), NMR; 2Z0D
# ATG4B (A) with rat LC3B (B); 9E8P ATG3 (A, B) with GABARAP (C, D)
DEFAULT = [
    ('pf_predicted_pf3d7_1209400_ciao2a', 'pf_reference_ciao1_ciao2a', None),
    ('pf_predicted_pf3d7_1356400_ppp2r1b', 'pf_reference_ppp2r3b_ppp2r1b', ('4i5l', 'A', 'B', 'B')),
    ('lb_predicted_atg8c_1_nbr1', 'lb_reference_map1lc3b_nbr1', ('2l8j', 'A', 'A', 'B')),
    ('lb_predicted_atg8c_1_atg4b', 'lb_reference_map1lc3b_atg4b', ('2z0d', 'A', 'B', 'B')),
    ('bm035_predicted_cmu_040530_snx15', None, None),
    ('bm012_predicted_cyth3_cnksr1', None, None),
    ('bm042_predicted_nas_14_3_rbm22', None, None),
    ('lb_predicted_atg8c_1_atg3', 'lb_reference_map1lc3b_atg3', ('9e8p', 'A', 'B', 'C+D')),
]


def model_file(job, model=0):
    job = job.lower().replace('-', '_')
    found = sorted(glob.glob(os.path.join(RESULTS, 'folds_*', job, f'*_model_{model}.cif')))
    if not found:
        raise CmdException(f'no model {model} of {job} in {RESULTS}')

    return found[-1]


def af_list(pattern=''):
    '''Prints the jobs with models whose names contain the pattern.'''
    jobs = sorted({os.path.basename(os.path.dirname(f)) for f in
                   glob.glob(os.path.join(RESULTS, 'folds_*', '*', '*_model_0.cif'))})
    for job in jobs:
        if pattern.lower() in job:
            print(job)


def object_name(job, prefix=''):
    '''
    `01_pred_pf3d7_1209400_ciao2a` for job pf_predicted_pf3d7_1209400_ciao2a with prefix
    01: the case study or benchmark link dropped, and the kind shortened.
    '''
    parts = job.lower().replace('-', '_').split('_')
    kind = next(i for i, part in enumerate(parts) if part in KINDS)
    name = '_'.join([KINDS[parts[kind]]] + parts[kind + 1:])

    return f'{prefix}_{name}' if prefix else name


def af_load(job, model=0, min_plddt=MIN_PLDDT, name=None):
    '''
    Loads a model of a job, its chains coloured, its interface as sticks and its residues
    below the pLDDT cutoff hidden.
    '''
    name = name or object_name(job)
    cmd.load(model_file(job, int(model)), name)
    cmd.hide('everything', name)
    cmd.show('cartoon', f'{name} and b > {float(min_plddt)}')
    cmd.color(PARASITE, f'{name} and chain A and elem C')
    cmd.color(HOST, f'{name} and chain B and elem C')
    interface = f'{name}_interface'
    a, b = f'({name} and chain A)', f'({name} and chain B)'
    cmd.select(interface, f'(byres ({a} within {INTERFACE} of {b})) or '
                          f'(byres ({b} within {INTERFACE} of {a}))')
    cmd.show('sticks', f'{interface} and b > {float(min_plddt)} and not name N+C+O')
    util.cnc(name)
    cmd.deselect()
    print(f'{name}: {cmd.count_atoms(interface + " and name CA")} interface residues')

    return name


def af_compare(predicted, reference, model=0, min_plddt=MIN_PLDDT, prefix=''):
    '''
    The predicted complex and its reference, superposed on their shared host partner where
    the reference binds it, so that the parasite protein and the host ortholog it stands
    in for can be compared.
    '''
    pred = af_load(predicted, model, min_plddt, object_name(predicted, prefix))
    ref = af_load(reference, model, min_plddt, object_name(reference, prefix))
    cmd.color(ORTHOLOG, f'{ref} and chain A and elem C')
    cmd.set('cartoon_transparency', 0.6, f'{ref} and chain B')
    cmd.hide('sticks', f'{ref} and chain B')
    # the host partner's confident residues near where the reference binds it, so that a
    # partner whose domains or tails lie differently in the two models is superposed at
    # the binding site; it is the same sequence in both models, so residues pair by number
    def confident(selection):
        space = {'found': set()}
        cmd.iterate(f'({selection}) and name CA and b > {CONFIDENT}', 'found.add(resv)',
                    space=space)
        return space['found']

    partner = confident(f'{pred} and chain B')
    site = partner & confident(f'{ref} and chain B within {SITE} of ({ref} and chain A)')
    on = f'its residues of pLDDT > {CONFIDENT} within {SITE} Å of the reference ortholog'
    if len(site) < 10:
        site = partner & confident(f'{ref} and chain B')
        on = f'its residues of pLDDT > {CONFIDENT}'
    resi = '+'.join(map(str, sorted(site)))
    rmsd = cmd.pair_fit(f'{ref} and chain B and name CA and resi {resi}',
                        f'{pred} and chain B and name CA and resi {resi}')
    print(f'host partner superposed on {on}, RMSD {rmsd:.2f} Å')
    # framed on what is shown, not the hidden low-pLDDT tails
    cmd.orient(f'({pred} or {ref}) and visible')

    return [pred, ref]


def af_pdb(entry, chain, target, others='', name=None):
    '''
    Fetches a PDB entry, keeps its first model and the chain superposed with the other
    chains given (of several, the one touching it most), and superposes that chain on the
    target, the superposed chain grey and the others violet, their contacts as sticks.
    '''
    entry = entry.lower()
    name = name or f'pdb_{entry}'
    os.makedirs(PDB, exist_ok=True)
    cmd.set('fetch_path', PDB)
    cmd.fetch(entry, f'_{name}', type='cif', async_=0)
    cmd.create(name, f'_{name} and state 1', 1, 1)
    cmd.delete(f'_{name}')
    others = [c for c in others.split('+') if c]
    if len(others) > 1:
        touching = {c: cmd.count_atoms(f'{name} and chain {c} within {INTERFACE} of '
                                       f'({name} and chain {chain})') for c in others}
        others = [max(touching, key=touching.get)]
    cmd.remove(f'{name} and not chain {"+".join([chain] + others)}')
    cmd.remove(f'{name} and (solvent or hetatm)')
    rmsd = cmd.align(f'{name} and chain {chain} and name CA', f'({target}) and name CA')[0]
    cmd.hide('everything', name)
    cmd.show('cartoon', name)
    cmd.color(HOST, f'{name} and chain {chain} and elem C')
    cmd.set('cartoon_transparency', 0.5, f'{name} and chain {chain}')
    if others:
        cmd.color(EXPERIMENTAL, f'{name} and chain {"+".join(others)} and elem C')
        cmd.show('sticks', f'byres ({name} and chain {"+".join(others)} within {INTERFACE} '
                           f'of ({name} and chain {chain})) and not name N+C+O')
    util.cnc(name)
    print(f'{name}: chain {chain} superposed on {target}, RMSD {rmsd:.2f} Å')

    return name


def af_plddt(selection='all'):
    '''Colours by pLDDT, stored as B-factors: >90, 70-90, 50-70, <50.'''
    cmd.color('0x0053d6', f'({selection}) and b > 90')
    cmd.color('0x65cbf3', f'({selection}) and b > 70 and b < 90.01')
    cmd.color('0xffdb13', f'({selection}) and b > 50 and b < 70.01')
    cmd.color('0xff7d45', f'({selection}) and b < 50.01')


cmd.extend('af_list', af_list)
cmd.extend('af_load', af_load)
cmd.extend('af_compare', af_compare)
cmd.extend('af_pdb', af_pdb)
cmd.extend('af_plddt', af_plddt)


def open_default():
    '''Loads the default set, one scene per complex: Page Up/Down steps through them.'''
    cmd.bg_color('white')
    cmd.set('cartoon_fancy_helices', 1)
    shown = {}
    for i, (predicted, reference, structure) in enumerate(DEFAULT, 1):
        prefix = f'{i:02d}'
        try:
            objects = (af_compare(predicted, reference, prefix=prefix) if reference else
                       [af_load(predicted, name=object_name(predicted, prefix))])
            if structure:
                entry, chain, target, others = structure
                objects.append(af_pdb(entry, chain, f'{objects[0]} and chain {target}',
                                      others, f'{prefix}_pdb_{entry}'))
            shown[prefix] = objects
        except CmdException as error:
            print(error)
    for prefix, objects in shown.items():
        cmd.disable('all')
        for name in objects:
            cmd.enable(name)
        cmd.orient('visible')
        cmd.scene(prefix, 'store')
    if shown:
        cmd.scene(next(iter(shown)), 'recall')


open_default()
