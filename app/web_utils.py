import base64
import os
import numpy as np
import pandas as pd
import streamlit as st
from streamlit_option_menu import option_menu
import utils

# menu entries as (label, bootstrap icon, path relative to the app root); the page file
# numbers no longer match this order
PAGES = [('Home', 'house', 'OrthoHPI_Home.py'),
         ('Parasites of a host', 'bug', 'pages/2_Compare_Parasites.py'),
         ('Multi-host parasites', 'person', 'pages/4_Hosts_of_a_Parasite.py'),
         ('Host-parasite network', 'share', 'pages/1_Predicted_Host-Parasite_PPIs.py'),
         ('About', 'info-circle', 'pages/3_About.py')]

TITLE = 'OrthoHPI 2.0'
SUBTITLE = 'Orthology Prediction of Host-Parasite PPIs'


def show_header(current=None, index=None):
    '''
    Draws the chrome every page opens with -- the name of the app, what it does, and the
    menu across the top -- and switches to whichever entry is picked, so that a page
    only has to say which one it is.
    '''
    st.markdown(f"<h1 style='text-align: center; color: #023858;'>{TITLE}</h1>",
                unsafe_allow_html=True)
    st.markdown(f"<h3 style='text-align: center; color: #2b8cbe;'>{SUBTITLE}</h3>",
                unsafe_allow_html=True)

    labels = [label for label, _, _ in PAGES]
    if index is None:
        index = labels.index(current) if current in labels else 0

    selected = option_menu(
        menu_title=None,
        options=labels,
        icons=[icon for _, icon, _ in PAGES],
        menu_icon="cast",
        default_index=index,
        orientation="horizontal",
    )

    for label, _, path in PAGES:
        if selected == label and selected != labels[index]:
            st.switch_page(path)

    return selected

# width figures are drawn at until the browser has reported the real page width
DEFAULT_PAGE_WIDTH = 1000
PAGE_WIDTH_SCRIPT = """(() => {
    const block = document.querySelector('[data-testid="stMainBlockContainer"]');
    if (!block) return null;
    const style = getComputedStyle(block);
    return block.clientWidth - parseFloat(style.paddingLeft) - parseFloat(style.paddingRight);
})()"""
# the gap streamlit leaves between two columns of a row, in pixels
COLUMN_GAP = 16


def page_width(default=DEFAULT_PAGE_WIDTH):
    '''Width of the page in pixels, on the screen the app is being read on.'''
    try:
        from streamlit_extras.eval_javascript import eval_javascript
        measured = int(eval_javascript(PAGE_WIDTH_SCRIPT, key='page_width'))
    except Exception:
        return default

    # a page narrower than a phone is one that has not finished laying itself out
    return measured if measured >= 320 else default


def column_width(columns, **kwargs):
    '''Width in pixels of one of `columns` equal columns of the page.'''
    return (page_width(**kwargs) - COLUMN_GAP * (columns - 1)) / columns


def get_data_dir():
    return st.session_state.get('data_dir', 'data')

def get_config_file():
    return st.session_state.get('config_file', 'config.yml')


def load_predictions(data_dir, config_file=None):
    '''
    Every predicted interaction that could take place: the ones whose host protein is
    expressed in a tissue the parasite is known to infect.
    '''
    return _load_predictions(data_dir, config_file or get_config_file())


@st.cache_data(show_spinner=False)
def _load_predictions(data_dir, config_file):
    predictions = utils.read_parquet_file(input_file=f'{data_dir}/predictions.parquet')
    predictions['weight'] = predictions['weight'].astype(float)

    return keep_infected_tissues(predictions, data_dir, config_file)


def keep_infected_tissues(predictions, data_dir, config_file):
    '''
    Drops the interactions whose host protein is not expressed anywhere the parasite is.
    '''
    tissue_file = os.path.join(data_dir, 'tissues_cell_types.parquet')
    if not os.path.exists(tissue_file):
        return predictions

    config = utils.read_config(config_file)
    tissues = utils.read_parquet_file(input_file=tissue_file)
    tissues = tissues.rename({'Gene': 'target'}, axis=1)[['target', 'Tissue']]
    expressed = tissues.drop_duplicates().groupby('Tissue')['target'].apply(frozenset)

    names = config['tissues']
    reachable = {}
    for taxid, parasite in config['parasites'].items():
        proteins = set()
        for tissue in parasite['tissues']:
            proteins |= expressed.get(names[tissue].lower(), frozenset())
        reachable[str(taxid)] = proteins

    keep = [target in reachable.get(str(taxid), ())
            for taxid, target in zip(predictions['taxid1'], predictions['target'])]

    return predictions[keep]


@st.cache_data(show_spinner=False)
def load_tissue_annotation(data_dir):
    '''
    The tissue and cell-type table, which pipeline/main.py writes for exactly the
    proteins that came through the secretome, tissue and DeepLoc filters.
    '''
    tissue_file = os.path.join(data_dir, 'tissues_cell_types.parquet')
    if not os.path.exists(tissue_file):
        return None

    return utils.read_parquet_file(input_file=tissue_file)


@st.cache_data(show_spinner=False)
def load_eligible_proteins(data_dir):
    '''
    The proteins the pipeline's filters passed, over every species, as pipeline/main.py
    writes them (or scripts/build_eligible_proteins.py for a data directory built before
    that file existed).
    '''
    eligible_file = os.path.join(data_dir, 'eligible_proteins.parquet')
    if not os.path.exists(eligible_file):
        return None

    return utils.read_parquet_file(input_file=eligible_file)


def load_proteome_sizes(data_dir):
    '''
    How many proteins STRING holds for every species, before any filter, as
    scripts/build_proteome_sizes.py writes them; None where the file was not built.
    '''
    sizes_file = os.path.join(data_dir, 'proteome_sizes.parquet')
    if not os.path.exists(sizes_file):
        return None

    return utils.read_parquet_file(input_file=sizes_file)


def filtered_pool(data_dir, taxids, niche=None):
    '''
    The proteins of these species that the pipeline had to work with: the ones that came
    through every filter, whether or not an interaction was predicted for them.
    '''
    eligible = load_eligible_proteins(data_dir)
    if eligible is not None:
        rows = eligible['taxid'].isin(taxids)
        column = niche_pool_column(niche)
        if column is not None and column in eligible.columns:
            rows = rows & eligible[column]
        return set(eligible.loc[rows, 'protein'])

    tissues = load_tissue_annotation(data_dir)
    if tissues is None:
        return set()

    return set(tissues.loc[tissues['Gene'].str.split('.').str[0].isin(taxids), 'Gene'])


def infected_tissue_proteins(data_dir, config, parasite_taxid):
    '''
    The host proteins annotated to a tissue this parasite is known to infect, over every
    host at once.
    '''
    tissues = load_tissue_annotation(data_dir)
    if tissues is None:
        return set()

    mapped = config['tissues']
    infected = {mapped[t].lower() for t in config['parasites'][int(parasite_taxid)]['tissues']}

    return set(tissues.loc[tissues['Tissue'].isin(infected), 'Gene'])


@st.cache_data(show_spinner=False)
def get_host_predictions(data_dir, host_taxids):
    '''
    Predictions against the selected host species. host_taxids is a one-taxid tuple of
    strings so the selector and its callers retain the same interface.
    '''
    predictions = load_predictions(data_dir)

    return predictions[predictions['taxid2'].isin(host_taxids)]


@st.cache_data(show_spinner=False)
def load_protein_annotations(data_dir):
    '''
    Descriptive protein names, keyed by STRING id, written by
    pipeline/build_protein_annotations.py.
    '''
    input_file = os.path.join(data_dir, 'protein_annotations.parquet')
    if not os.path.exists(input_file):
        return {}
    annotations = utils.read_parquet_file(input_file=input_file)

    return dict(annotations[['protein', 'description']].values)


@st.cache_data(show_spinner=False)
def load_deeploc_localisations(data_dir):
    '''
    Where DeepLoc 2 predicts each protein of the predictions sits, keyed by STRING id,
    written by pipeline/build_deeploc_localisations.py.
    '''
    input_file = os.path.join(data_dir, 'deeploc_localisations.parquet')
    if not os.path.exists(input_file):
        return pd.DataFrame(columns=['protein', 'localizations', 'signals',
                                     'membrane_types'] + list(DEEPLOC_SCORES.values()))

    return utils.read_parquet_file(input_file=input_file)


CELL_MEMBRANE = 'Cell membrane'
EXTRACELLULAR = 'Extracellular'
CYTOPLASM = 'Cytoplasm'
NUCLEUS = 'Nucleus'
# BOTH_SURFACE where only the surface pair is read, SEVERAL where all four classes are
BOTH_SURFACE = 'Both'
SEVERAL = 'Several'
NOT_SURFACE = 'Neither'

# DeepLoc 2 Accurate-model thresholds, as the pipeline filters with. DeepLoc2/deeploc2.py
# label_threshold is read at i+1, so Extracellular is labels[2]; see docs/deeploc.md
DEEPLOC_CUTOFFS = {EXTRACELLULAR: 0.61728516, CELL_MEMBRANE: 0.56464844,
                   CYTOPLASM: 0.47612305, NUCLEUS: 0.50136719}
DEEPLOC_SCORES = {EXTRACELLULAR: 'extracellular', CELL_MEMBRANE: 'cell_membrane',
                  CYTOPLASM: 'cytoplasm', NUCLEUS: 'nucleus'}
# a host protein is kept for any of the four classes, a parasite protein by the secretome
# filter on the surface pair alone
SURFACE_CLASSES = (EXTRACELLULAR, CELL_MEMBRANE)
HOST_CLASSES = (EXTRACELLULAR, CELL_MEMBRANE, CYTOPLASM, NUCLEUS)

# two shades of one hue per pair: blue for the surface classes, orange for the intracellular
# ones; purple for several classes, grey for none
LOCALISATION_COLORS = {EXTRACELLULAR: '#a6bddb', CELL_MEMBRANE: '#045a8d',
                       CYTOPLASM: '#fdae6b', NUCLEUS: '#a63603',
                       BOTH_SURFACE: '#756bb1', SEVERAL: '#756bb1',
                       NOT_SURFACE: '#d9d9d9'}


def niche_classes(niche):
    '''
    The classes a host protein of this parasite is read on: the four of HOST_CLASSES for
    a parasite with an intracellular stage, the surface pair for one without.
    '''
    return HOST_CLASSES if niche == 'Intracellular' else SURFACE_CLASSES


def classify_localisation(localisations, classes=HOST_CLASSES, multiple=SEVERAL):
    '''
    Which localization class DeepLoc puts each protein in, on the same cut-offs the
    pipeline filtered it with.
    '''
    read = [c for c in classes if DEEPLOC_SCORES[c] in localisations]
    called = {c: localisations[DEEPLOC_SCORES[c]] > DEEPLOC_CUTOFFS[c] for c in read}
    if not called:
        return pd.Series(NOT_SURFACE, index=localisations.index)

    # np.select takes the first condition that holds, so several classes is tested before
    # any single one
    over_several = sum(called.values()) > 1

    return pd.Series(np.select([over_several] + [called[c] for c in read],
                               [multiple] + list(read), default=NOT_SURFACE),
                     index=localisations.index)


# where a parasite sits relative to the host cell, i.e. which host proteins it can reach;
# independent of `multicellular`
UNKNOWN_NICHE = 'Unknown'
NICHE_ORDER = ['Extracellular', 'Intracellular']
# greys, so the niche is not read as a third categorical split beside the clades and the
# DeepLoc classes
NICHE_COLORS = {'Extracellular': '#c7c7c7', 'Intracellular': '#3d3d3d',
                UNKNOWN_NICHE: '#f0f0f0'}
NICHE_TITLE = 'intracellular / extracellular'


def parasite_niche(config, taxid):
    '''
    The niche config.yml records for one parasite, keyed by taxid rather than by label,
    so a page holding predictions can look it up from the column they carry.
    '''
    parasite = config.get('parasites', {}).get(int(taxid), {})

    return str(parasite.get('niche', '')).strip().lower() or None


def niche_pool_column(niche):
    '''
    The eligible_proteins column holding the half of the pool a niche reaches, as
    pipeline/main.py names it, or None where there is no niche to narrow by.
    '''
    if not niche or str(niche).capitalize() not in NICHE_ORDER:
        return None

    return f'reachable_{str(niche).strip().lower()}'


def get_niches(config):
    '''
    The niche of every parasite of the configuration, keyed by the label the predictions
    name it with, so a figure can annotate a column without reading the config itself.
    '''
    niches = {}
    for parasite in config.get('parasites', {}).values():
        niche = str(parasite.get('niche', '')).capitalize()
        niches[parasite['label']] = niche if niche in NICHE_ORDER else UNKNOWN_NICHE

    return niches


def get_host_groups(config, predictions, include_rodents=False):
    '''Maps each host species offered in the app to its taxid.'''
    predicted = set(predictions['taxid2'])
    groups = {}
    for taxid, host in config['hosts'].items():
        taxid = str(taxid)
        if taxid in predicted:
            groups[host['label']] = [taxid]

    rodent_taxids = ['10116', '10090']
    if include_rodents and all(taxid in predicted for taxid in rodent_taxids):
        groups['Rodent'] = rodent_taxids

    return {group: groups[group] for group in sorted(groups)}


# Streamlit drops the state of a widget the current page does not draw, so the host is kept
# in a key of our own that every page seeds from and writes back to
HOST_STATE_KEY = 'selected_host'
HOST_WIDGET_KEY = 'selected_host_widget'
NO_HOST = '<select>'


def host_selector(config, predictions, label='Select a host', include_rodents=False):
    '''Draws the host selectbox, shared across pages through HOST_STATE_KEY.'''
    groups = get_host_groups(config, predictions, include_rodents=include_rodents)
    options = [NO_HOST] + list(groups)

    current = st.session_state.get(HOST_STATE_KEY, NO_HOST)
    if current not in options:
        # the stored host has no predictions in this data_dir (snapshot entrypoints), so
        # fall back
        current = NO_HOST

    # seeded through the key rather than `index`: a changing `index` makes Streamlit treat
    # it as a new widget and drop the click that changed the host
    if st.session_state.get(HOST_WIDGET_KEY) not in options:
        st.session_state[HOST_WIDGET_KEY] = current

    selected = st.selectbox(label, options, key=HOST_WIDGET_KEY)
    st.session_state[HOST_STATE_KEY] = selected

    return selected, tuple(groups.get(selected, ()))


def filter_tissues(config, df):
    source = df['taxid1'].unique()[0]
    mapped_tissues = config['tissues']
    tissues = [mapped_tissues[t].lower() for t in config['parasites'][int(source)]['tissues']]
    df = df[df['Tissue'].isin(tissues)]
    
    return df


# minimum atlas expression for a host protein to count in a cell type
CELL_TYPE_NTPM_CUTOFF = 1.0


def keep_expressed_cell_types(df, cutoff=CELL_TYPE_NTPM_CUTOFF):
    '''
    The rows of a tissue-annotated frame where a host protein has expression above the
    absolute cell-type cutoff.
    '''
    return df[df['nTPM'] > cutoff]


def count_ticks(figure, largest, axis='x', **kwargs):
    '''Put whole numbers on an axis that counts things.'''
    update = figure.update_xaxes if axis == 'x' else figure.update_yaxes
    if largest <= 10:
        update(dtick=1, tickformat='d', **kwargs)
    else:
        update(tickformat='d', **kwargs)


# (name, link, logo file) in pipeline order; scripts/build_footer_logos.py puts the files on
# one background and height
LOGOS = [('EggNOG', 'http://eggnog6.embl.de/', 'eggnog.png'),
         ('STRING', 'https://string-db.org/', 'string.png'),
         ('TISSUES', 'https://tissues.jensenlab.org/', 'tissues.png'),
         ('Human Protein Atlas', 'https://www.proteinatlas.org/', 'hpa.png'),
         ('Pig Cell Atlas', 'https://dreamapp.biomed.au.dk/pigatlas/',
          'pca.png'),
         ('Tabula Muris Senis',
          'https://cellxgene.cziscience.com/collections/0b9d8a04-bb9d-44da-aa27-705bb65b54eb',
          'tabula_muris_senis.png'),
         ('Gene Ontology', 'https://geneontology.org/', 'go.png'),
         ('AlphaFold', 'https://deepmind.google/science/alphafold/', 'deepmind.png'),
         ('AlphaFold Protein Structure Database', 'https://alphafold.ebi.ac.uk/',
          'embl-ebi.svg')]

LOGO_DIR = os.path.join('images', 'logos')

# the files are written at twice this so they stay sharp on a retina screen
LOGO_HEIGHT = 40

REPOSITORY = 'https://github.com/Multiomics-Analytics-Group/OrthoHPI2.0'

MIME_TYPES = {'.png': 'image/png', '.svg': 'image/svg+xml'}


@st.cache_data(show_spinner=False)
def logo_source(filename):
    '''
    A logo as a data uri, so that it can be given a link. st.image cannot be wrapped in
    one, and a plain <img src="images/..."> is not served by Streamlit at all.
    '''
    with open(os.path.join(LOGO_DIR, filename), 'rb') as f:
        encoded = base64.b64encode(f.read()).decode()

    mime = MIME_TYPES[os.path.splitext(filename)[1]]

    return f'data:{mime};base64,{encoded}'


def footer():
    st.write("Developed with data from:")

    logos = ''.join(
        f'<a href="{link}" target="_blank" title="{name}">'
        f'<img src="{logo_source(filename)}" alt="{name}" height="{LOGO_HEIGHT}"></a>'
        for name, link, filename in LOGOS)

    st.markdown(
        '<div style="display: flex; flex-wrap: wrap; align-items: center; gap: 1.5rem 2.5rem;'
        f' margin-bottom: 1rem;">{logos}</div>', unsafe_allow_html=True)

    st.caption('Tissue names follow the BRENDA Tissue Ontology, and the localisation of '
               'the host proteins is predicted with DeepLoc 2.')

    st.markdown(f'Code available at: [{REPOSITORY.split("//")[1]}]({REPOSITORY})')