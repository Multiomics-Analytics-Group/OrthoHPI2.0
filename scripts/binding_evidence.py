'''
Human protein interactions in BioGRID and IntAct by what kind of binding they show:
direct (two-hybrid, reconstituted complex, co-crystal structure and other assays of two
purified or fused proteins, or an enzymatic reaction of one on the other), co-complex
(affinity capture, co-fractionation, co-purification), or proximity only (proximity
labelling, and in IntAct also resonance energy transfer). BioGRID's kind is read from
its experimental system, IntAct's from its interaction type. Downloads both databases
into data/downloads/ the first time; used by scripts/build_alphafold_direct_jobs.py.
'''
import datetime
import gzip
import os
import shutil
import urllib.request
import zipfile

import pandas as pd

BIOGRID_URL = ('https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/'
               'BIOGRID-ALL-LATEST.tab3.zip')
INTACT_URL = 'https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip'
HUMAN = 9606
# most direct first
KINDS = ['direct', 'co-complex', 'proximity']
# BioGRID's physical experimental systems by what they show
BIOGRID_KINDS = {
    'direct': {'Two-hybrid', 'Reconstituted Complex', 'Co-crystal Structure', 'Far Western',
               'FRET', 'PCA', 'Protein-peptide', 'Biochemical Activity'},
    'co-complex': {'Affinity Capture-MS', 'Affinity Capture-Western',
                   'Affinity Capture-Luminescence', 'Affinity Capture-RNA',
                   'Co-fractionation', 'Co-purification'},
    'proximity': {'Proximity Label-MS'},
}


def download(url, path):
    if not os.path.exists(path):
        print(f'downloading {url}')
        # BioGRID refuses Python's default user agent
        request = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(request) as response, open(f'{path}.part', 'wb') as f:
            shutil.copyfileobj(response, f)
        os.rename(f'{path}.part', path)

    return path


def pairs_of(a, b):
    return [tuple(sorted(p)) for p in zip(a.str.upper(), b.str.upper())]


def biogrid(downloads):
    '''Human physical interactions in BioGRID, with their kind, and the release.'''
    path = download(BIOGRID_URL, os.path.join(downloads, 'BIOGRID-ALL-LATEST.tab3.zip'))
    with zipfile.ZipFile(path) as archive:
        name = archive.namelist()[0]
        with archive.open(name) as f:
            rows = pd.read_csv(f, sep='\t', low_memory=False, usecols=[
                'Official Symbol Interactor A', 'Official Symbol Interactor B',
                'Experimental System', 'Experimental System Type',
                'Organism ID Interactor A', 'Organism ID Interactor B', 'Publication Source'])
    rows = rows[(rows['Experimental System Type'] == 'physical') &
                (rows['Organism ID Interactor A'] == HUMAN) &
                (rows['Organism ID Interactor B'] == HUMAN)]
    kinds = {system: kind for kind, systems in BIOGRID_KINDS.items() for system in systems}
    rows = pd.DataFrame({
        'pair': pairs_of(rows['Official Symbol Interactor A'], rows['Official Symbol Interactor B']),
        'method': rows['Experimental System'].to_numpy(),
        'kind': rows['Experimental System'].map(kinds).to_numpy(),
        'publication': rows['Publication Source'].str.lower().to_numpy()})

    return rows.assign(source='BioGRID'), name.removesuffix('.tab3.txt')


def intact_kind(kind):
    '''IntAct's interaction type as a kind; None for colocalization and the like.'''
    if kind == 'proximity':
        return 'proximity'
    if kind in ('direct interaction', 'protein cleavage', 'covalent binding',
                'disulfide bond') or kind.endswith('reaction'):
        return 'direct'
    if kind in ('physical association', 'association'):
        return 'co-complex'

    return None


def intact(downloads):
    '''
    Human interactions in IntAct, with their kind, and the date of the download. The
    human rows are kept in data/downloads/intact_human.parquet, as reading the full file
    takes minutes. IntAct names no release, so it goes by the date of its file.
    '''
    path = download(INTACT_URL, os.path.join(downloads, 'intact.zip'))
    with zipfile.ZipFile(path) as archive:
        release = 'IntAct ' + datetime.date(*archive.getinfo('intact.txt').date_time[:3]).isoformat()
    cache = os.path.join(downloads, 'intact_human.parquet')
    if os.path.exists(cache) and os.path.getmtime(cache) > os.path.getmtime(path):
        rows = pd.read_parquet(cache)
        return rows.assign(kind=rows['type'].map(intact_kind)), release

    print('reading IntAct')
    name = lambda column: column.str.extract(r'\(([^()]*)\)\s*$', expand=False).fillna('')
    gene = lambda column: column.str.extract(r'uniprotkb:([^|()]+)\(gene name\)', expand=False)
    chunks = []
    with zipfile.ZipFile(path) as archive:
        with archive.open('intact.txt') as f:
            for chunk in pd.read_csv(f, sep='\t', usecols=[4, 5, 6, 8, 9, 10, 11], header=0,
                                     names=['aliases A', 'aliases B', 'method', 'publications',
                                            'taxid A', 'taxid B', 'type'],
                                     dtype=str, chunksize=500_000, quoting=3):
                human = (chunk['taxid A'].str.startswith(f'taxid:{HUMAN}(', na=False) &
                         chunk['taxid B'].str.startswith(f'taxid:{HUMAN}(', na=False))
                chunk = chunk[human]
                a, b = gene(chunk['aliases A']), gene(chunk['aliases B'])
                keep = a.notna() & b.notna()
                chunk, a, b = chunk[keep], a[keep], b[keep]
                method, kind = name(chunk['method']), name(chunk['type'])
                chunks.append(pd.DataFrame({
                    'pair': pairs_of(a, b), 'method': method.to_numpy(),
                    'type': kind.to_numpy(),
                    'publication': chunk['publications'].str.extract(
                        r'(pubmed:\d+)', expand=False).to_numpy()}))
    rows = pd.concat(chunks, ignore_index=True)
    rows['pair'] = rows['pair'].map('|'.join)
    rows.to_parquet(cache)

    return rows.assign(kind=rows['type'].map(intact_kind)), release


def human_members(downloads, groups, symbols):
    '''The gene symbols of each group's human members, from the EggNOG members file.'''
    members = {}
    with gzip.open(os.path.join(downloads, '2759_members.tsv.gz'), 'rt') as f:
        for line in f:
            fields = line.split('\t')
            if fields[1] in groups:
                members[fields[1]] = {symbols[p].upper() for p in fields[4].split(',')
                                      if p.startswith(f'{HUMAN}.') and p in symbols}

    return members


def interactions(downloads):
    '''
    Human interactions of known kind in BioGRID and IntAct, pairs of upper-case symbols
    joined by |, and the two releases.
    '''
    from_biogrid, biogrid_release = biogrid(downloads)
    from_intact, intact_release = intact(downloads)
    from_biogrid['pair'] = from_biogrid['pair'].map('|'.join)
    rows = pd.concat([from_biogrid, from_intact.assign(source='IntAct')], ignore_index=True)

    return rows[rows['kind'].notna()], f'{biogrid_release}; {intact_release}'


def most_direct(kinds):
    return next((k for k in KINDS if k in kinds), 'none found')
