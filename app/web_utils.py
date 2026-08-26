import base64
import os
import numpy as np
import pandas as pd
import streamlit as st
from streamlit_option_menu import option_menu
import utils

# The pages of the app in the order they are offered, each as (label, bootstrap icon,
# path to switch to). The app reads at four levels -- every host, one host, one parasite
# across the hosts it infects, one host and one parasite -- and the menu follows that
# order. Each icon shows what the page gives back rather than what it is asked about, so
# the two directional pages read as a pair and stay apart at menu size. The paths are
# relative to the app root, which is what st.switch_page takes; the file numbers of the
# pages no longer match this order (the network page keeps the URL it was published
# under) and nothing reads them, since the Streamlit sidebar is hidden in css/style.css.
PAGES = [('Home', 'house', 'OrthoHPI_Home.py'),
         ('Parasites of a host', 'bug', 'pages/2_Compare_Parasites.py'),
         ('Hosts of a parasite', 'person', 'pages/4_Hosts_of_a_Parasite.py'),
         ('Host-parasite network', 'share', 'pages/1_Predicted_Host-Parasite_PPIs.py'),
         ('About', 'info-circle', 'pages/3_About.py')]

# the name of the app and what it does, drawn above the menu on every page
TITLE = 'OrthoHPI 2.0'
SUBTITLE = 'Orthology Prediction of Host-Parasite PPIs'


def show_header(current=None, index=None):
    '''
    Draws the chrome every page opens with -- the name of the app, what it does, and the
    menu across the top -- and switches to whichever entry is picked, so that a page only
    has to say which one it is.

    The title goes above the menu and is drawn here rather than by each page: the pages
    carried a copy of the same two lines each, below their own menu, which put the
    navigation above the name of the thing being navigated.

    :param str current: label of the page drawing the menu; the entry it highlights and
                        the one click that is not a navigation
    :param int index: position of that page in PAGES, for the snapshot entrypoints that
                      predate the labels
    :return: the label selected, which is `current` unless the page is being left
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

def get_data_dir():
    return st.session_state.get('data_dir', 'data')

def get_config_file():
    return st.session_state.get('config_file', 'config.yml')


def load_predictions(data_dir, config_file=None):
    '''
    Every predicted interaction that could take place: the ones whose host protein is
    expressed in a tissue the parasite is known to infect. Shared by the pages rather than
    reloaded on each of them, so the restriction is applied once and every page counts the
    same interactions.

    The config file is resolved here rather than inside the cache so that it is part of the
    cache key: a snapshot entrypoint pointing at another configuration infects other
    tissues and keeps other interactions.

    :param str data_dir: directory holding predictions.parquet
    :param str config_file: configuration to read the infected tissues from, defaulting to
                            the one the session was started with
    :return: predictions dataframe
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

    The pipeline filters the host proteins by tissue already, but on the union of the
    tissues every parasite infects, so a gut parasite keeps a protein expressed only in
    brain. Applying it per parasite is what makes an interaction on any page a claim about
    something that could happen: without it a parasite of one tissue carries predictions
    against proteins it never meets -- for the narrowest of them, 99 of every 100.

    A data directory with no tissue table is left alone rather than emptied, which is what
    the snapshot directories need.

    :param predictions: predictions dataframe
    :param str data_dir: directory holding tissues_cell_types.parquet
    :param str config_file: configuration naming the tissues each parasite infects
    :return: the predictions whose host protein is expressed where its parasite is
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
    The tissue and cell-type table, which pipeline/main.py writes for exactly the proteins
    that came through the secretome, tissue and DeepLoc filters. Host proteins only: the
    parasites are filtered on their secretome and never enter it.

    :param str data_dir: directory holding tissues_cell_types.parquet
    :return: the table, or None where the directory does not carry one
    '''
    tissue_file = os.path.join(data_dir, 'tissues_cell_types.parquet')
    if not os.path.exists(tissue_file):
        return None

    return utils.read_parquet_file(input_file=tissue_file)


@st.cache_data(show_spinner=False)
def load_eligible_proteins(data_dir):
    '''
    The proteins the pipeline's filters passed, over every species, as
    pipeline/main.py writes them (or scripts/build_eligible_proteins.py for a data
    directory built before that file existed).

    :param str data_dir: directory holding eligible_proteins.parquet
    :return: the table, or None where the directory does not carry one
    '''
    eligible_file = os.path.join(data_dir, 'eligible_proteins.parquet')
    if not os.path.exists(eligible_file):
        return None

    return utils.read_parquet_file(input_file=eligible_file)


def filtered_pool(data_dir, taxids):
    '''
    The proteins of these species that the pipeline had to work with: the ones that came
    through every filter, whether or not an interaction was predicted for them.

    This is the set a network is drawn from, so it is the background any enrichment of
    that network has to be read against. Tested against the whole proteome instead, a
    network reports the filters that built it -- its host proteins were kept for being
    surface or extracellular and its parasite proteins for being secreted, and those are
    the processes that come out enriched whichever parasite is asked about.

    A data directory written before eligible_proteins.parquet existed falls back to the
    tissue table, which covers the hosts alone; its parasites are then left on the
    proteome, which is what they were tested against before.

    :param str data_dir: directory holding eligible_proteins.parquet
    :param tuple taxids: taxids as strings, of either side
    :return: set of STRING protein ids, empty where the directory carries neither table
    '''
    eligible = load_eligible_proteins(data_dir)
    if eligible is not None:
        return set(eligible.loc[eligible['taxid'].isin(taxids), 'protein'])

    tissues = load_tissue_annotation(data_dir)
    if tissues is None:
        return set()

    return set(tissues.loc[tissues['Gene'].str.split('.').str[0].isin(taxids), 'Gene'])


def infected_tissue_proteins(data_dir, config, parasite_taxid):
    '''
    The host proteins annotated to a tissue this parasite is known to infect, over every
    host at once.

    The pipeline's tissue filter is not parasite-specific -- it keeps a host protein
    expressed in a tissue *any* parasite of the config infects -- so a host protein can
    carry a predicted interaction with a parasite that never reaches the tissue it was
    kept for. This is the set that says which ones do not.

    :param dict config: parsed configuration, naming the tissues each parasite infects
    :param parasite_taxid: taxid of the parasite, as anything int() takes
    :return: set of STRING protein ids
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
    Predictions against one host group. host_taxids is a tuple of taxids (as strings)
    so grouped hosts -- rat and mouse under Rodent -- are pooled into a single view;
    the predictions keep their own species label and taxid.

    :param str data_dir: directory holding predictions.parquet
    :param tuple host_taxids: taxids of the selected host group, as strings
    :return: predictions dataframe of that host group
    '''
    predictions = load_predictions(data_dir)

    return predictions[predictions['taxid2'].isin(host_taxids)]


@st.cache_data(show_spinner=False)
def load_protein_annotations(data_dir):
    '''
    Descriptive protein names, keyed by STRING id, written by
    pipeline/build_protein_annotations.py. The snapshot data directories predate that
    script, so a missing file only leaves the protein descriptions empty.

    :param str data_dir: directory holding protein_annotations.parquet
    :return: {STRING id: descriptive protein name}
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
    written by pipeline/build_deeploc_localisations.py. The pipeline reads the same
    predictions as a filter and keeps only the verdict, so this is what is left of the
    probabilities behind it. The snapshot data directories predate that script, so a
    missing file only leaves the localisations out of the figures that show them.

    :param str data_dir: directory holding deeploc_localisations.parquet
    :return: dataframe of protein, localizations, signals, membrane_types,
             extracellular, cell_membrane; empty if the file is not there
    '''
    input_file = os.path.join(data_dir, 'deeploc_localisations.parquet')
    if not os.path.exists(input_file):
        return pd.DataFrame(columns=['protein', 'localizations', 'signals', 'membrane_types',
                                     'extracellular', 'cell_membrane'])

    return utils.read_parquet_file(input_file=input_file)


# the DeepLoc 2 (Accurate) probabilities a host protein is kept as surface-exposed on,
# the same cut-offs the pipeline filters with (pipeline/main.py). They are read in the app
# to say which of the two reasons a protein was kept for, not to filter anything again
DEEPLOC_MEMBRANE_CUTOFF = 0.56464844
DEEPLOC_EXTRACELLULAR_CUTOFF = 0.61728516
# the two ways of being surface-exposed, and what a protein over neither cut-off is called
CELL_MEMBRANE = 'Cell membrane'
EXTRACELLULAR = 'Extracellular'
NOT_SURFACE = 'Neither'


def classify_surface(localisations):
    '''
    Which of the two surface classes DeepLoc puts each protein in.

    A protein over both cut-offs is called extracellular: the few there are (two of the
    human ones, the laminins LAMA1 and LAMB2, each barely over the membrane cut-off) do
    not earn a third class in every figure that shows one. A protein over neither is
    NOT_SURFACE -- on the host side that is a host the filter was never run for, or a
    data directory built with other cut-offs; on the parasite side it is ordinary, since
    parasite proteins are selected by the secretome filter and not by DeepLoc.

    :param localisations: DeepLoc table, as load_deeploc_localisations returns it
    :return: series of class names, aligned to the rows of the table
    '''
    return pd.Series(np.select([localisations['extracellular'] >= DEEPLOC_EXTRACELLULAR_CUTOFF,
                                localisations['cell_membrane'] >= DEEPLOC_MEMBRANE_CUTOFF],
                               [EXTRACELLULAR, CELL_MEMBRANE], default=NOT_SURFACE),
                     index=localisations.index)


def get_host_groups(config, predictions):
    '''
    Maps each host group offered in the app to the taxids it covers. Hosts sharing a
    `group` in the config (rat and mouse -> Rodent) are presented as one option, hosts
    without one stand alone under their label. Groups with no predicted interaction are
    dropped so the selector never offers an empty host.

    :param dict config: parsed configuration
    :param predictions: predictions dataframe, used to drop groups without predictions
    :return: {group label: [taxid as str, ...]} sorted by group label
    '''
    predicted = set(predictions['taxid2'])
    groups = {}
    for taxid, host in config['hosts'].items():
        taxid = str(taxid)
        if taxid in predicted:
            groups.setdefault(host.get('group', host['label']), []).append(taxid)

    return {group: groups[group] for group in sorted(groups)}


# The host chosen on any page applies to all of them. Streamlit drops the state of a
# widget that the current page does not draw, so the choice is kept in a key of our
# own (HOST_STATE_KEY) that nothing else touches, and the selectbox on each page is
# seeded from it and writes back.
HOST_STATE_KEY = 'selected_host'
HOST_WIDGET_KEY = 'selected_host_widget'
NO_HOST = '<select>'


def host_selector(config, predictions, label='Select a host'):
    '''
    Draws the host selectbox, shared across pages through HOST_STATE_KEY.

    :param dict config: parsed configuration
    :param predictions: predictions dataframe, used to build the options
    :param str label: label shown above the selectbox
    :return: (group label, tuple of taxids as str); NO_HOST and () when nothing is chosen
    '''
    groups = get_host_groups(config, predictions)
    options = [NO_HOST] + list(groups)

    current = st.session_state.get(HOST_STATE_KEY, NO_HOST)
    if current not in options:
        # the stored host has no predictions in this data_dir (the snapshot entrypoints
        # each point at their own), so fall back rather than raise
        current = NO_HOST

    # The widget is seeded through its own key rather than through `index`. Streamlit
    # identifies a keyless widget by its arguments, so a changing `index` makes it a
    # different widget: it dropped the click that changed the host and re-rendered the
    # previous one, and the host only changed on the second click. Seeding is skipped
    # once the widget holds a valid choice, so the user's selection is not overwritten.
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

# The sources the predictions are built from, as (name, link, logo file), in the order the
# pipeline reaches them: the orthologous groups, the interactions transferred along them,
# then what the host protein has to be to keep a prediction, and the structures shown of
# one. The files are prepared by scripts/build_footer_logos.py, which puts them all on the
# same background and the same height -- st.columns of a fixed width left them at anything
# between 19 and 65 pixels tall, since the wordmarks are nothing like the same shape.
LOGOS = [('EggNOG', 'http://eggnog6.embl.de/', 'eggnog.png'),
         ('STRING', 'https://string-db.org/', 'string.png'),
         ('TISSUES', 'https://tissues.jensenlab.org/', 'tissues.png'),
         ('Human Protein Atlas', 'https://www.proteinatlas.org/', 'hpa.png'),
         ('Gene Ontology', 'https://geneontology.org/', 'go.png'),
         ('AlphaFold', 'https://deepmind.google/science/alphafold/', 'deepmind.png'),
         ('AlphaFold Protein Structure Database', 'https://alphafold.ebi.ac.uk/',
          'embl-ebi.svg')]

LOGO_DIR = os.path.join('images', 'logos')

# what the footer draws them at; scripts/build_footer_logos.py writes the files at twice
# this, so they stay sharp on a retina screen
LOGO_HEIGHT = 40

REPOSITORY = 'https://github.com/Multiomics-Analytics-Group/OrthoHPI2.0'

MIME_TYPES = {'.png': 'image/png', '.svg': 'image/svg+xml'}


@st.cache_data(show_spinner=False)
def logo_source(filename):
    '''
    A logo as a data uri, so that it can be given a link. st.image cannot be wrapped in
    one, and a plain <img src="images/..."> is not served by Streamlit at all.

    :param str filename: name of the file in images/logos
    :return: the file as a data uri
    '''
    with open(os.path.join(LOGO_DIR, filename), 'rb') as f:
        encoded = base64.b64encode(f.read()).decode()

    mime = MIME_TYPES[os.path.splitext(filename)[1]]

    return f'data:{mime};base64,{encoded}'


def footer():
    st.write("Developed with data from:")

    # a wrapping row rather than fixed columns, which squeezed the logos into each other
    # on a narrow window instead of moving them onto a second line
    logos = ''.join(
        f'<a href="{link}" target="_blank" title="{name}">'
        f'<img src="{logo_source(filename)}" alt="{name}" height="{LOGO_HEIGHT}"></a>'
        for name, link, filename in LOGOS)

    st.markdown(
        '<div style="display: flex; flex-wrap: wrap; align-items: center; gap: 1.5rem 2.5rem;'
        f' margin-bottom: 1rem;">{logos}</div>', unsafe_allow_html=True)

    # neither has a logo to put in the row above, and leaving them out credited the tissue
    # vocabulary and the localisation filter to nobody
    st.caption('Tissue names follow the BRENDA Tissue Ontology, and the localisation of '
               'the host proteins is predicted with DeepLoc 2.')

    st.markdown(f'Code available at: [{REPOSITORY.split("//")[1]}]({REPOSITORY})')