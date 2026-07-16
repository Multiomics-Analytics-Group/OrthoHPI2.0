import os
import urllib.request
from stmol import makeobj


def get_alphafold_structure(query_proteins={}):
    structures = {}
    url = 'https://alphafold.ebi.ac.uk/files/AF-query_protein-F1-model_v4.pdb'
    web_url = 'https://alphafold.ebi.ac.uk/entry/query_protein'
    for query_protein in query_proteins:
        pdb_filename = None
        uniprot_id = query_proteins[query_protein]
        try:
            if uniprot_id is not None:
                request = urllib.request.Request(url.replace('query_protein', uniprot_id))
                pdb_filename = os.path.join('data/tmp', query_protein+'_output_structure.pdb')
                if not os.path.isfile(pdb_filename):
                    with open(pdb_filename, 'w') as out:
                        with urllib.request.urlopen(request) as response:
                            res = response.read().decode('utf-8')
                            out.write(res)
                structures[query_protein] = (pdb_filename,
                                        url.replace('query_protein', uniprot_id),
                                        web_url.replace('query_protein', uniprot_id))
        except Exception as e:
            print(e)

    return structures

def generate_mol_structure(pdb_file):
    with open(pdb_file) as ifile:
        content = ifile.read()

    xyzview = makeobj(content, molformat='pdb', style='cartoon', background='black')
    xyzview.setStyle({'cartoon':{'color':'spectrum'}})
    
    return xyzview