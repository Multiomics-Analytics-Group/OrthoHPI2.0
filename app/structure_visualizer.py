import os
import json
import urllib.request
from functools import lru_cache
import py3Dmol


api_url = 'https://alphafold.ebi.ac.uk/api/prediction/query_protein'
web_url = 'https://alphafold.ebi.ac.uk/entry/query_protein'


@lru_cache(maxsize=None)
def get_alphafold_url(uniprot_id):
    '''
    Returns the url of the current AlphaFold model file for a UniProt id. The model
    files are versioned and older versions are removed (model_v4 is gone), so the
    url is resolved through the API instead of being built from a fixed version.

    :param str uniprot_id: UniProt accession to look up
    :return: url of the pdb file for the latest model of that protein
    '''
    request = urllib.request.Request(api_url.replace('query_protein', uniprot_id))
    with urllib.request.urlopen(request, timeout=30) as response:
        predictions = json.loads(response.read().decode('utf-8'))

    if not predictions:
        raise ValueError(f'no AlphaFold prediction for {uniprot_id}')

    return predictions[0]['pdbUrl']


def download_structure(url, pdb_filename):
    '''
    Downloads a structure to pdb_filename. The file is written to a temporary path
    first so a failed download does not leave an empty file behind, which would
    later be mistaken for a cached structure.
    '''
    request = urllib.request.Request(url)
    with urllib.request.urlopen(request, timeout=60) as response:
        content = response.read().decode('utf-8')

    tmp_filename = pdb_filename + '.part'
    with open(tmp_filename, 'w') as out:
        out.write(content)

    os.replace(tmp_filename, pdb_filename)


def get_alphafold_structure(query_proteins={}, output_dir='data/tmp'):
    '''
    Downloads the AlphaFold structure of each query protein.

    :param dict query_proteins: protein name --> UniProt id
    :param str output_dir: directory where the downloaded structures are cached
    :return: dict protein name --> (pdb file, pdb url, AlphaFold entry url), where
             the file and url are None when no structure could be retrieved
    '''
    structures = {}
    os.makedirs(output_dir, exist_ok=True)
    for query_protein in query_proteins:
        uniprot_id = query_proteins[query_protein]
        pdb_filename = None
        pdb_url = None
        entry_url = web_url.replace('query_protein', str(uniprot_id))
        try:
            if uniprot_id is not None:
                pdb_url = get_alphafold_url(uniprot_id)
                pdb_filename = os.path.join(output_dir, query_protein+'_output_structure.pdb')
                if not os.path.isfile(pdb_filename) or os.path.getsize(pdb_filename) == 0:
                    download_structure(pdb_url, pdb_filename)
        except Exception as e:
            print(f'No AlphaFold structure for {query_protein} ({uniprot_id}): {e}')
            pdb_filename = None
            pdb_url = None

        structures[query_protein] = (pdb_filename, pdb_url, entry_url)

    return structures

def generate_mol_structure(pdb_file, height=460):
    '''
    Builds the 3Dmol viewer of a structure. The viewer is built to fill the width it is
    given rather than the fixed 640px py3Dmol defaults to, so that it is not cut off when
    embedded in a column narrower than that. The view is also zoomed out a little from
    the fit 3Dmol computes, which otherwise leaves elongated proteins touching the edges
    of the canvas.

    :param str pdb_file: path of the pdb file to show
    :param int height: height of the viewer in pixels
    :return: py3Dmol view of the structure
    '''
    with open(pdb_file) as ifile:
        content = ifile.read()

    xyzview = py3Dmol.view(width='100%', height=height)
    xyzview.addModel(content, 'pdb')
    xyzview.setStyle({'cartoon': {'color': 'spectrum'}})
    xyzview.setBackgroundColor('black')
    xyzview.zoomTo()
    xyzview.zoom(0.85)

    return xyzview
