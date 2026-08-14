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
    :return: dict protein name --> (pdb file, pdb url, AlphaFold entry url, reason),
             where the file and url are None when no structure could be retrieved and
             reason says why, and reason is None when one was
    '''
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
            # 404 is the database answering that it holds no model for this protein,
            # which about one parasite protein in five hits and is not an error. Anything
            # else went wrong on the way and should not be reported as a missing model.
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

# AlphaFold writes the per-residue confidence of its model (pLDDT) into the B-factor
# column of the pdb file, so colouring the cartoon by that column shows which parts of
# the model are to be trusted. The colours are the ones the AlphaFold database itself
# uses, from very low to very high confidence, and the range is cut at the two outer
# thresholds of its scale: 3Dmol blends between the four rather than banding them, so
# the bands read as a gradient from orange through yellow and light blue to dark blue.
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

    :return: html of the legend, to be written with st.markdown(unsafe_allow_html=True)
    '''
    swatches = ''.join(
        f'<span style="display:inline-block;width:0.75em;height:0.75em;'
        f'background:{color};border-radius:2px;margin:0 0.3em 0 1em;"></span>{label}'
        for color, label in PLDDT_BANDS)

    return (f'<div style="font-size:0.8em;color:#555555;">Model confidence{swatches}</div>')


def generate_mol_structure(pdb_file, height=460):
    '''
    Builds the 3Dmol viewer of a structure, with the cartoon coloured by the per-residue
    confidence of the model. The viewer is built to fill the width it is given rather
    than the fixed 640px py3Dmol defaults to, so that it is not cut off when embedded in
    a column narrower than that. The view is also zoomed out a little from the fit 3Dmol
    computes, which otherwise leaves elongated proteins touching the edges of the canvas.

    :param str pdb_file: path of the pdb file to show
    :param int height: height of the viewer in pixels
    :return: py3Dmol view of the structure
    '''
    with open(pdb_file) as ifile:
        content = ifile.read()

    xyzview = py3Dmol.view(width='100%', height=height)
    xyzview.addModel(content, 'pdb')
    xyzview.setStyle({'cartoon': {'colorscheme': {'prop': 'b', 'gradient': 'linear',
                                                  'colors': PLDDT_COLORS,
                                                  'min': PLDDT_MIN, 'max': PLDDT_MAX}}})
    # white rather than the black py3Dmol starts from: the confident parts of the model
    # are the dark blue end of the scale, which does not stand out against black
    xyzview.setBackgroundColor('white')
    xyzview.zoomTo()
    xyzview.zoom(0.85)

    return xyzview
