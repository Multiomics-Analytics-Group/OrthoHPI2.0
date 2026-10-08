'''
What kind of experiment supports each link of the AlphaFold benchmark
(scripts/build_alphafold_benchmark_jobs.py), from the human members of its two
orthologous groups in BioGRID and IntAct: direct (two-hybrid, reconstituted complex,
co-crystal structure and other assays of two purified or fused proteins, or an
enzymatic reaction of one on the other), co-complex (affinity capture,
co-fractionation, co-purification), or proximity only (proximity labelling, and in
IntAct also resonance energy transfer). BioGRID's kind is read from its experimental
system, IntAct's from its interaction type. A link is labelled by
the most direct kind either database shows for any human pair of its groups. Downloads
both databases into data/downloads/ the first time. Writes
paper/tables/alphafold_benchmark_evidence.csv.
'''
import argparse
import datetime
import gzip
import os
import re
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


def build(benchmark, data_dir, out):
    downloads = os.path.join(data_dir, 'downloads')
    links = pd.read_csv(benchmark)
    links = links[links['kind'] == 'predicted']
    eligible = pd.read_parquet(os.path.join(data_dir, 'eligible_proteins.parquet'))
    human = eligible[eligible['taxid'].astype(str) == str(HUMAN)]
    symbols = dict(zip(human['protein'], human['name']))
    members = human_members(downloads, set(links['group A']) | set(links['group B']), symbols)
    known, releases = interactions(downloads)
    by_pair = {pair: rows for pair, rows in known.groupby('pair')}

    rows = []
    for link in links.to_dict('records'):
        a, b = members.get(link['group A'], set()), members.get(link['group B'], set())
        pairs = {'|'.join(sorted((x, y))) for x in a for y in b if x != y}
        found = pd.concat([by_pair[p] for p in pairs if p in by_pair] or [known.iloc[:0]])
        row = {'job': link['job'], 'link': link['link'], 'evidence': link['evidence'],
               'score': link['score'], 'human A': len(a), 'human B': len(b),
               'pairs with evidence': found['pair'].nunique(),
               'publications': found['publication'].dropna().nunique(),
               'kind': most_direct(set(found['kind']))}
        for source in ('BioGRID', 'IntAct'):
            hits = found[found['source'] == source]
            row[f'kind {source}'] = most_direct(set(hits['kind']))
            row[f'methods {source}'] = '; '.join(sorted(set(hits['method'])))
        rows.append(row)

    table = pd.DataFrame(rows).assign(releases=releases)
    table.to_csv(out, index=False)
    print(f'{releases}; {len(table)} links, '
          f'{(table["human A"] * table["human B"] == 0).sum()} with a group lacking human members')
    print(pd.crosstab(table['kind'], table['evidence'], margins=True))
    print(pd.crosstab(table['kind BioGRID'], table['kind IntAct'], margins=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--benchmark', default='paper/alphafold/benchmark.csv')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--out', default='paper/tables/alphafold_benchmark_evidence.csv')
    args = parser.parse_args()
    build(args.benchmark, args.data_dir, args.out)
