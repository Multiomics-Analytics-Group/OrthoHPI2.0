import web_utils
import streamlit as st

config = web_utils.load_config(web_utils.get_config_file())
hosts = ', '.join(host['label'].split(' (')[1].rstrip(')') for host in config['hosts'].values())

st.title('About')

st.markdown(f'''
OrthoHPI predicts protein–protein interactions between {len(config['parasites'])} parasites
and their hosts ({hosts}). Interactions known between orthology groups in STRING are
transferred to the parasite and host proteins in those groups, and kept only where the
parasite protein can reach the host protein.

#### How predictions are made

1. **Parasite proteins**: proteins predicted by DeepLoc 2 to be secreted. Unicellular
   parasites also keep their cell-membrane proteins.
2. **Host proteins**: proteins expressed in a tissue the parasite infects (TISSUES database) and located where the parasite can reach them (DeepLoc 2): the cell surface and
   extracellular space, plus the cytoplasm and nucleus for intracellular parasites.
3. **Interactions**: a parasite and a host protein are paired when their EggNOG 6 orthology
   groups interact in STRING v12.0 with experimental or database evidence ≥ 0.7.

Cell-type expression from the Human Protein Atlas, the Pig Cell Atlas and Tabula Muris
Senis annotates the predictions. Rat has no cell-type data.

#### Contact

Ida K.S. Meitil (idamei@dtu.dk) and Alberto Santos (albsad@dtu.dk), Technical University of Denmark.
''')

with st.container():
    web_utils.footer()
