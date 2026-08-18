import os
import numpy as np
import pandas as pd
import streamlit as st
from streamlit_option_menu import option_menu
import utils

# The pages of the app in the order they are offered, each as (label, bootstrap icon,
# path to switch to). The app reads at four levels -- every host, one host, one parasite
# across the hosts it infects, one host and one parasite -- and the menu follows that
# order. The paths are relative to
# the app root, which is what st.switch_page takes; the file numbers of the pages no
# longer match this order (the network page keeps the URL it was published under) and
# nothing reads them, since the Streamlit sidebar is hidden in css/style.css.
PAGES = [('Home', 'house', 'OrthoHPI_Home.py'),
         ('Parasites in a host', 'diagram-2', 'pages/2_Compare_Parasites.py'),
         ('One parasite, several hosts', 'signpost-split', 'pages/4_Hosts_of_a_Parasite.py'),
         ('Predicted Host-parasite PPIs', 'diagram-3', 'pages/1_Predicted_Host-Parasite_PPIs.py'),
         ('About', 'chat-text', 'pages/3_About.py')]


def show_pages_menu(current=None, index=None):
    '''
    Draws the menu across the top of every page and switches to whichever entry is
    picked, so that a page only has to say which one it is.

    :param str current: label of the page drawing the menu; the entry it highlights and
                        the one click that is not a navigation
    :param int index: position of that page in PAGES, for the snapshot entrypoints that
                      predate the labels
    :return: the label selected, which is `current` unless the page is being left
    '''
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


@st.cache_data(show_spinner=False)
def load_predictions(data_dir):
    '''
    Every predicted interaction. Shared by the three pages rather than reloaded on each
    of them: the cache is keyed by the data directory, so the parquet is read once
    however the app is navigated.

    :param str data_dir: directory holding predictions.parquet
    :return: predictions dataframe
    '''
    predictions = utils.read_parquet_file(input_file=f'{data_dir}/predictions.parquet')
    predictions['weight'] = predictions['weight'].astype(float)

    return predictions


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

def footer():
    st.write("Developed with data from:")

    cols = st.columns(5)
    with cols[0]:
        st.image('images/eggnog.png', width=200)
    with cols[1]:
        st.image('images/string.png', width=200)
    with cols[2]:
        st.image('images/hpa.png', width=200)
    with cols[3]:
        st.image('images/tissues.png', width=200)
    with cols[4]:
        st.image('images/ebi.png', width=200)

    st.write("Code available at: https://github.com/Multiomics-Analytics-Group/OrthoHPI2.0")