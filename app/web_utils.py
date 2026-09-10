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
         ('Multi-host parasites', 'person', 'pages/4_Hosts_of_a_Parasite.py'),
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

# the width the figures that carry a size of their own are drawn for until the browser has
# said how wide the page really is. On the narrow side of the screens the app is read on, so
# that the first draw of a session is one that fits its column rather than one that is cut
# off at the edge of it: a figure too small for a moment is the better way to be wrong
DEFAULT_PAGE_WIDTH = 1000
# what the page is measured with, run in the browser: the inner width of the block the app
# is drawn in, which is the width the figures share between them, rather than the width of
# the window, which also holds the padding either side of that block
PAGE_WIDTH_SCRIPT = """(() => {
    const block = document.querySelector('[data-testid="stMainBlockContainer"]');
    if (!block) return null;
    const style = getComputedStyle(block);
    return block.clientWidth - parseFloat(style.paddingLeft) - parseFloat(style.paddingRight);
})()"""
# the gap streamlit leaves between two columns of a row, in pixels
COLUMN_GAP = 16


def page_width(default=DEFAULT_PAGE_WIDTH):
    '''
    Width of the page in pixels, on the screen the app is being read on.

    Streamlit runs on the server and is told nothing about the window, so a figure that has
    to know how much room it has -- one whose cells must come out square, say -- has to ask
    the browser, which answers on a later run of the script. Nothing should depend on the
    answer arriving: a figure is drawn at `default` until it does, and redrawn at the
    measured width on the run it arrives on.

    :return: the width the browser reported, or `default` before it has reported one, in a
             browser that refused to answer, or outside a browser altogether (a test run)
    '''
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


def filtered_pool(data_dir, taxids, niche=None):
    '''
    The proteins of these species that the pipeline had to work with: the ones that came
    through every filter, whether or not an interaction was predicted for them.

    The host half of that pool depends on the niche of the parasite being asked about --
    the cytosolic and nuclear proteins were only ever open to the intracellular parasites
    -- so a niche narrows it to the column pipeline/main.py wrote for that niche. A pool
    built before those columns existed is returned whole, as it was read before.

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
    :param str niche: niche of the parasite the pool is a background for, as config.yml
                      records it; None leaves the pool at its union over the niches
    :return: set of STRING protein ids, empty where the directory carries neither table
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
    Predictions against the selected host species. host_taxids is a one-taxid tuple of
    strings so the selector and its callers retain the same interface.

    :param str data_dir: directory holding predictions.parquet
    :param tuple host_taxids: taxid of the selected host species, as a string
    :return: predictions dataframe of that host species
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
    :return: dataframe of protein, localizations, signals, membrane_types and the
             probability of each class of DEEPLOC_SCORES; empty if the file is not there
    '''
    input_file = os.path.join(data_dir, 'deeploc_localisations.parquet')
    if not os.path.exists(input_file):
        return pd.DataFrame(columns=['protein', 'localizations', 'signals',
                                     'membrane_types'] + list(DEEPLOC_SCORES.values()))

    return utils.read_parquet_file(input_file=input_file)


# the localization classes a parasite can reach: the two surface ones either kind of
# parasite meets, and the two an intracellular parasite reaches inside the host cell
CELL_MEMBRANE = 'Cell membrane'
EXTRACELLULAR = 'Extracellular'
CYTOPLASM = 'Cytoplasm'
NUCLEUS = 'Nucleus'
# the class of a protein called for more than one of them -- BOTH_SURFACE where only the
# two surface classes are being read, SEVERAL where all four are -- and the fallback for a
# protein called for none of the classes it is read on
BOTH_SURFACE = 'Both'
SEVERAL = 'Several'
NOT_SURFACE = 'Neither'

# DeepLoc 2's own per-class thresholds for the Accurate (ProtT5) model, the same values the
# pipeline filters both sides with -- the host proteins in pipeline/main.py and the parasite
# proteins in deeploc/build_secretome_fastas.py. Taken from DeepLoc2/deeploc2.py
# label_threshold, which carries one entry more than there are classes and is read at i+1,
# so Extracellular is labels[2] -> 0.61728516 and Cell membrane is labels[3] -> 0.56464844.
# docs/deeploc.md has the derivation.
DEEPLOC_CUTOFFS = {EXTRACELLULAR: 0.61728516, CELL_MEMBRANE: 0.56464844,
                   CYTOPLASM: 0.47612305, NUCLEUS: 0.50136719}
# the column of deeploc_localisations.parquet each class is scored on
DEEPLOC_SCORES = {EXTRACELLULAR: 'extracellular', CELL_MEMBRANE: 'cell_membrane',
                  CYTOPLASM: 'cytoplasm', NUCLEUS: 'nucleus'}
# The classes each side of an interaction is read on, mirroring the two filters that built
# it (pipeline/main.py DEEPLOC_NICHE_CLASSES, deeploc/build_secretome_fastas.py). A host
# protein is kept for any of four, the cytosol and the nucleus being open to a parasite
# with an intracellular stage; a parasite protein is kept by the secretome filter, which
# reads the surface pair alone, so its own cytosolic proteins say nothing about what it
# reaches its host with and are not a class it is split by.
SURFACE_CLASSES = (EXTRACELLULAR, CELL_MEMBRANE)
HOST_CLASSES = (EXTRACELLULAR, CELL_MEMBRANE, CYTOPLASM, NUCLEUS)

# The colour of each class, wherever a figure of the app shows one -- the bars of the home
# page and the strip beside the rows of the shared-interactor dot plot -- so that a class
# learned on one page is read on the next without a second key. Each pair is two shades of
# one hue: the blue the pages are headed in for the two surface classes, an orange for the
# two inside the cell, so a figure reads as two places a parasite can meet its host rather
# than as four unrelated categories, and the shades within a pair say it is one whole split
# up. A protein called for more than one class is neither shade of either and gets a purple
# of its own; one called for none of them, which is a species the filter never ran on, a
# grey. The dark shade of each pair is darker than the taxonomic group nearest it -- the
# strips carry Okabe-Ito, whose vermillion #D55E00 is the nearest thing to these oranges --
# so a bar is not read as the clade of the parasite it stands over.
LOCALISATION_COLORS = {EXTRACELLULAR: '#a6bddb', CELL_MEMBRANE: '#045a8d',
                       CYTOPLASM: '#fdae6b', NUCLEUS: '#a63603',
                       BOTH_SURFACE: '#756bb1', SEVERAL: '#756bb1',
                       NOT_SURFACE: '#d9d9d9'}


def niche_classes(niche):
    '''
    The classes a host protein of this parasite is read on: the four of HOST_CLASSES for a
    parasite with an intracellular stage, the surface pair for one without.

    A host protein can be called for a class its parasite's niche never gave it -- plenty
    of cell membrane proteins are cytosolic as well -- and reading it on that class would
    say the parasite meets it there, which is what apply_deeploc_filter did not allow. So
    the classes follow the niche, and a parasite whose niche is unknown is read on the
    surface pair, the narrower of the two, exactly as the filter treated it.

    :param str niche: niche of the parasite, as get_niches names it
    :return: the classes to read its host proteins on
    '''
    return HOST_CLASSES if niche == 'Intracellular' else SURFACE_CLASSES


def classify_localisation(localisations, classes=HOST_CLASSES, multiple=SEVERAL):
    '''
    Which localization class DeepLoc puts each protein in, on the same cut-offs the pipeline
    filtered it with. A protein called for more than one of them is `multiple`.

    The classes are the caller's, because the two sides of an interaction were filtered on
    different ones: HOST_CLASSES for a host protein, which an intracellular parasite reaches
    in the cytosol and the nucleus as well as at the surface, and SURFACE_CLASSES for a
    parasite protein, which the secretome filter kept for the surface pair alone. A class
    whose probability the table does not carry is dropped, so a snapshot data directory
    written before pipeline/build_deeploc_localisations.py kept the cytosolic
    probabilities is read on the classes it does carry rather than failing.

    NOT_SURFACE, called for none of them, is a protein of a species the filter was never
    run for, or one of the snapshot directories built with other cut-offs. It does not
    occur in a data directory this pipeline built and read on the classes that built it:
    every protein of the predictions is over the cut-off of at least one of them.

    Read from the probabilities and not from the `localizations` column beside them. That
    column names a class even for a protein that crosses no threshold at all -- DeepLoc
    falls back to whichever class came closest, by the largest probability minus its own
    threshold (DeepLoc2/deeploc2.py) -- so it would put proteins in a class the filters
    rejected and draw them below the cut-off line of every figure that shows one.

    :param localisations: DeepLoc table, as load_deeploc_localisations returns it
    :param classes: the classes to read the table on, HOST_CLASSES or SURFACE_CLASSES
    :param str multiple: what to call a protein over the cut-off of more than one of them:
                         BOTH_SURFACE where there are two classes, SEVERAL where there are
                         more, the word being what the legend of a figure carries
    :return: series of class names, aligned to the rows of the table
    '''
    read = [c for c in classes if DEEPLOC_SCORES[c] in localisations]
    called = {c: localisations[DEEPLOC_SCORES[c]] > DEEPLOC_CUTOFFS[c] for c in read}
    if not called:
        return pd.Series(NOT_SURFACE, index=localisations.index)

    # np.select takes the first condition that holds, so more than one class is tested for
    # before any single class, which would otherwise claim the protein
    over_several = sum(called.values()) > 1

    return pd.Series(np.select([over_several] + [called[c] for c in read],
                               [multiple] + list(read), default=NOT_SURFACE),
                     index=localisations.index)


# Where a parasite sits relative to the host cell, as `niche` records it in config.yml:
# which host proteins it is in a position to reach at all. It is not the same statement as
# `multicellular`, which is about the other side of the interface -- which of the parasite's
# own proteins are exposed to the host -- and the two are independent: Trichinella is
# multicellular and intracellular, Trypanosoma brucei unicellular and extracellular.
UNKNOWN_NICHE = 'Unknown'
# outside the host cell, through to inside it, which is the order the bands and the legend
# draw the values in
NICHE_ORDER = ['Extracellular', 'Intracellular']
# Greys, and deliberately not a hue of their own. The clades already have the Okabe-Ito set
# and the DeepLoc classes the blues, and a third categorical palette beside those two would
# be read as a third thing the bars are split into rather than as a fact about the parasite
# under them. A ramp from pale to dark also carries the order of the values, which a set of
# separate hues would not.
NICHE_COLORS = {'Extracellular': '#c7c7c7', 'Intracellular': '#3d3d3d',
                UNKNOWN_NICHE: '#f0f0f0'}
# what the band and its legend are called. `niche` is the config key, but it is not the
# term the literature uses for this split -- that is the two words themselves -- so the
# figures name the values rather than the key
NICHE_TITLE = 'intracellular / extracellular'


def parasite_niche(config, taxid):
    '''
    The niche config.yml records for one parasite, keyed by taxid rather than by label, so
    a page holding predictions can look it up from the column they carry.

    :param dict config: parsed configuration
    :param taxid: parasite taxid, as an int or the string the predictions carry
    :return: the niche as the config spells it, or None where it records none
    '''
    parasite = config.get('parasites', {}).get(int(taxid), {})

    return str(parasite.get('niche', '')).strip().lower() or None


def niche_pool_column(niche):
    '''
    The eligible_proteins column holding the half of the pool a niche reaches, as
    pipeline/main.py names it, or None where there is no niche to narrow by.

    :param str niche: niche as config.yml records it, or a display value of NICHE_ORDER
    :return: column name, or None
    '''
    if not niche or str(niche).capitalize() not in NICHE_ORDER:
        return None

    return f'reachable_{str(niche).strip().lower()}'


def get_niches(config):
    '''
    The niche of every parasite of the configuration, keyed by the label the predictions
    name it with, so a figure can annotate a column without reading the config itself.

    :param dict config: parsed configuration
    :return: {parasite label: one of NICHE_ORDER, or UNKNOWN_NICHE}
    '''
    niches = {}
    for parasite in config.get('parasites', {}).values():
        niche = str(parasite.get('niche', '')).capitalize()
        niches[parasite['label']] = niche if niche in NICHE_ORDER else UNKNOWN_NICHE

    return niches


def get_host_groups(config, predictions, include_rodents=False):
    '''
    Maps each host species offered in the app to its taxid. Hosts with no predicted
    interaction are dropped so the selector never offers an empty host.

    :param dict config: parsed configuration of individual host species
    :param predictions: predictions dataframe, used to drop groups without predictions
    :param bool include_rodents: add a combined rat/mouse option
    :return: {host species label: [taxid as str]} sorted by label
    '''
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


# The host chosen on any page applies to all of them. Streamlit drops the state of a
# widget that the current page does not draw, so the choice is kept in a key of our
# own (HOST_STATE_KEY) that nothing else touches, and the selectbox on each page is
# seeded from it and writes back.
HOST_STATE_KEY = 'selected_host'
HOST_WIDGET_KEY = 'selected_host_widget'
NO_HOST = '<select>'


def host_selector(config, predictions, label='Select a host', include_rodents=False):
    '''
    Draws the host selectbox, shared across pages through HOST_STATE_KEY.

    :param dict config: parsed configuration
    :param predictions: predictions dataframe, used to build the options
    :param str label: label shown above the selectbox
    :param bool include_rodents: offer a combined rat/mouse option
    :return: (group label, tuple of taxids as str); NO_HOST and () when nothing is chosen
    '''
    groups = get_host_groups(config, predictions, include_rodents=include_rodents)
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


# Minimum atlas expression for a host protein to count in a cell type.
CELL_TYPE_NTPM_CUTOFF = 1.0


def keep_expressed_cell_types(df, cutoff=CELL_TYPE_NTPM_CUTOFF):
    '''
    The rows of a tissue-annotated frame where a host protein has expression above the
    absolute cell-type cutoff.

    HPA rows are retained in the generated annotation from nTPM > 0. Applying this
    cutoff when counting or filtering by cell type avoids treating negligible expression as
    cell-type presence.

    Rows the HPA gives no cell type -- every host but human, its single cell data being
    human only -- carry no nTPM either and drop out.

    :param dataframe df: rows carrying target, Tissue and nTPM
    :param float cutoff: minimum nTPM a cell type needs, exclusive
    :return: the subset of the rows
    '''
    return df[df['nTPM'] > cutoff]


def count_ticks(figure, largest, axis='x', **kwargs):
    '''
    Put whole numbers on an axis that counts things.

    Left to itself plotly divides a short range into halves, so a figure whose tallest bar
    is three proteins is ruled at 0.5, 1.0, 1.5 -- and half a protein is not a quantity.
    Above ten its own steps are whole numbers already and are left alone, since one tick
    per unit would crowd the axis.

    :param figure: the figure to rule
    :param largest: the largest count drawn
    :param str axis: 'x' or 'y'
    :param kwargs: passed on to the axis, so a caller can rule it in one call
    '''
    update = figure.update_xaxes if axis == 'x' else figure.update_yaxes
    if largest <= 10:
        update(dtick=1, tickformat='d', **kwargs)
    else:
        update(tickformat='d', **kwargs)


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