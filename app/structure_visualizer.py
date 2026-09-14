import os
import json
import urllib.request
import urllib.error
from functools import lru_cache
import py3Dmol


api_url = 'https://alphafold.ebi.ac.uk/api/prediction/query_protein'
web_url = 'https://alphafold.ebi.ac.uk/entry/query_protein'


@lru_cache(maxsize=None)
def get_alphafold_url(uniprot_id):
    '''Returns the url of the current AlphaFold model file for a UniProt id.'''
    request = urllib.request.Request(api_url.replace('query_protein', uniprot_id))
    with urllib.request.urlopen(request, timeout=30) as response:
        predictions = json.loads(response.read().decode('utf-8'))

    if not predictions:
        raise ValueError(f'no AlphaFold prediction for {uniprot_id}')

    return predictions[0]['pdbUrl']


def download_structure(url, pdb_filename):
    '''Downloads a structure to pdb_filename.'''
    request = urllib.request.Request(url)
    with urllib.request.urlopen(request, timeout=60) as response:
        content = response.read().decode('utf-8')

    tmp_filename = pdb_filename + '.part'
    with open(tmp_filename, 'w') as out:
        out.write(content)

    os.replace(tmp_filename, pdb_filename)


def get_alphafold_structure(query_proteins={}, output_dir='data/tmp'):
    '''Downloads the AlphaFold structure of each query protein.'''
    structures = {}
    os.makedirs(output_dir, exist_ok=True)
    for query_protein in query_proteins:
        uniprot_id = query_proteins[query_protein]
        pdb_filename = None
        pdb_url = None
        entry_url = web_url.replace('query_protein', str(uniprot_id))
        reason = None
        try:
            if uniprot_id is None:
                reason = 'This protein has no UniProt accession to look a model up by'
            else:
                pdb_url = get_alphafold_url(uniprot_id)
                pdb_filename = os.path.join(output_dir, query_protein+'_output_structure.pdb')
                if not os.path.isfile(pdb_filename) or os.path.getsize(pdb_filename) == 0:
                    download_structure(pdb_url, pdb_filename)
        except urllib.error.HTTPError as e:
            # 404 means no model for this protein, which about one parasite protein in five
            # hits
            if e.code == 404:
                reason = 'The AlphaFold database holds no model for this protein'
            else:
                reason = f'The AlphaFold database could not be read ({e})'
            print(f'No AlphaFold structure for {query_protein} ({uniprot_id}): {e}')
            pdb_filename = None
            pdb_url = None
        except Exception as e:
            reason = f'The AlphaFold database could not be reached ({e})'
            print(f'No AlphaFold structure for {query_protein} ({uniprot_id}): {e}')
            pdb_filename = None
            pdb_url = None

        structures[query_protein] = (pdb_filename, pdb_url, entry_url, reason)

    return structures

# AlphaFold writes pLDDT into the B-factor column; colours are the AlphaFold database's own
PLDDT_BANDS = [('#FF7D45', 'Very low (pLDDT < 50)'),
               ('#FFDB13', 'Low (50 - 70)'),
               ('#65CBF3', 'Confident (70 - 90)'),
               ('#0053D6', 'Very high (> 90)')]
PLDDT_COLORS = [color for color, _ in PLDDT_BANDS]
PLDDT_MIN = 50
PLDDT_MAX = 90


def plddt_legend():
    '''
    The key to the colours generate_mol_structure paints the model with, which mean
    nothing without one.
    '''
    swatches = ''.join(
        f'<span style="display:inline-block;width:0.75em;height:0.75em;'
        f'background:{color};border-radius:2px;margin:0 0.3em 0 1em;"></span>{label}'
        for color, label in PLDDT_BANDS)

    return (f'<div style="font-size:0.8em;color:#555555;">Model confidence{swatches}</div>')


def generate_mol_structure(pdb_file, height=460):
    '''
    Builds the 3Dmol viewer of a structure, with the cartoon coloured by the per-residue
    confidence of the model.
    '''
    with open(pdb_file) as ifile:
        content = ifile.read()

    xyzview = py3Dmol.view(width='100%', height=height)
    xyzview.addModel(content, 'pdb')
    xyzview.setStyle({'cartoon': {'colorscheme': {'prop': 'b', 'gradient': 'linear',
                                                  'colors': PLDDT_COLORS,
                                                  'min': PLDDT_MIN, 'max': PLDDT_MAX}}})
    # the confident parts of the model are dark blue, which does not stand out against black
    xyzview.setBackgroundColor('white')
    xyzview.zoomTo()
    xyzview.zoom(0.85)

    return xyzview
