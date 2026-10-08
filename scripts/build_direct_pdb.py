'''
Which human references of the direct-interaction benchmark
(scripts/build_alphafold_direct_jobs.py) have an experimental structure of the two
proteins together in the PDB, and whether one was released before AlphaFold 3's
training cutoff, so that the model may have learnt the complex. The proteins are mapped
to reviewed UniProt entries by gene name, and the PDB searched for entries holding both.
Writes paper/tables/alphafold_direct_pdb.csv.
'''
import argparse
import json
import urllib.parse
import urllib.request

import pandas as pd

# the last PDB release date in AlphaFold 3's training set (Abramson et al. 2024)
TRAINING_CUTOFF = '2021-09-30'


def get(url, body=None):
    request = urllib.request.Request(url, data=json.dumps(body).encode() if body else None,
                                     headers={'Content-Type': 'application/json'})
    with urllib.request.urlopen(request, timeout=60) as response:
        text = response.read()

    return json.loads(text) if text else None


def uniprot(symbol):
    '''The reviewed human UniProt accession of a gene name.'''
    query = urllib.parse.quote(f'gene_exact:{symbol} AND organism_id:9606 AND reviewed:true')
    found = get(f'https://rest.uniprot.org/uniprotkb/search?query={query}&fields=accession&format=json')
    accessions = [r['primaryAccession'] for r in found['results']]

    return accessions[0] if accessions else None


def entries(a, b):
    '''PDB entries holding both accessions, with their method and release date.'''
    def holds(accession):
        return {'type': 'terminal', 'service': 'text', 'parameters': {
            'attribute': 'rcsb_polymer_entity_container_identifiers.'
                         'reference_sequence_identifiers.database_accession',
            'operator': 'exact_match', 'value': accession}}

    found = get('https://search.rcsb.org/rcsbsearch/v2/query', {
        'query': {'type': 'group', 'logical_operator': 'and', 'nodes': [holds(a), holds(b)]},
        'return_type': 'entry', 'request_options': {'return_all_hits': True}})
    ids = [hit['identifier'] for hit in found['result_set']] if found else []
    if not ids:
        return []
    details = get('https://data.rcsb.org/graphql', {'query': '{entries(entry_ids: %s) '
                  '{rcsb_id exptl {method} rcsb_accession_info {initial_release_date}}}'
                  % json.dumps(ids)})

    return [(e['rcsb_id'], e['exptl'][0]['method'],
             e['rcsb_accession_info']['initial_release_date'][:10])
            for e in details['data']['entries']]


def build(jobs_file, out):
    references = pd.read_csv(jobs_file)
    references = references[references['kind'] == 'reference']
    symbols = sorted(set(references['name A']) | set(references['name B']))
    accessions = {s: uniprot(s) for s in symbols}
    rows = []
    for ref in references.to_dict('records'):
        a, b = accessions[ref['name A']], accessions[ref['name B']]
        found = entries(a, b) if a and b else []
        dates = sorted(date for _, _, date in found)
        rows.append({'job': ref['job'], 'link': ref['link'], 'name A': ref['name A'],
                     'name B': ref['name B'], 'uniprot A': a, 'uniprot B': b,
                     'structures': len(found),
                     'entries': ' '.join(f'{i} ({m.lower()}, {d})' for i, m, d in found),
                     'first released': dates[0] if dates else None,
                     'before training cutoff': bool(dates) and dates[0] <= TRAINING_CUTOFF})
    table = pd.DataFrame(rows)
    table.to_csv(out, index=False)
    print(f"{len(table)} references; {(table['structures'] > 0).sum()} with a structure, "
          f"{table['before training cutoff'].sum()} released by {TRAINING_CUTOFF}; "
          f"unmapped: {sorted(s for s, a in accessions.items() if a is None)}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--jobs', default='paper/alphafold/direct.csv')
    parser.add_argument('--out', default='paper/tables/alphafold_direct_pdb.csv')
    args = parser.parse_args()
    build(args.jobs, args.out)
