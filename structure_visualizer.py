import os
from stmol import makeobj
import urllib.request

def get_alphafold_structure(query_proteins={}):
    structures = {}
    url = 'https://alphafold.ebi.ac.uk/files/AF-query_protein-F1-model_v4.pdb'

    for query_protein in query_proteins:    
        uniprot_id = query_proteins[query_protein]
        try:
            request = urllib.request.Request(url.replace('query_protein', uniprot_id))
            pdb_filename = os.path.join('data/tmp', query_protein+'_output_structure.pdb')
            if not os.path.isfile(pdb_filename):
                with open(pdb_filename, 'w') as out:
                    with urllib.request.urlopen(request) as response:
                        res = response.read().decode('utf-8')
                        out.write(res)
        except Exception as e:
            pdb_filename = None
            print(e)

        structures[query_protein] = (pdb_filename, url.replace('query_protein', uniprot_id))

    return structures

def generate_mol_structure(pdb_file):
    with open(pdb_file) as ifile:
        content = ifile.read()

    xyzview = makeobj(content, molformat='pdb', style='cartoon', background='white')
    xyzview.setStyle({'cartoon':{'color':'spectrum'}})
    
    return xyzview