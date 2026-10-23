import base64
import os
import numpy as np
import pandas as pd
import streamlit as st
import textwrap
import utils

# menu entries as (label, material icon, path relative to the app root); the page file
# numbers no longer match this order
PAGES = [('Home', ':material/home:', 'views/Home.py'),
         ('Host-parasite network', ':material/share:', 'views/1_Predicted_Host-Parasite_PPIs.py'),
         ('Parasites of a host', ':material/bug_report:', 'views/2_Compare_Parasites.py'),
         ('Multi-host parasites', ':material/person:', 'views/4_Hosts_of_a_Parasite.py'),
         ('About', ':material/info:', 'views/3_About.py')]

# pixels one line of a wrapped axis name takes, at plotly's default font
LABEL_LINE = 15


def wrap_label(label, width, lines=4):
    '''
    An axis name broken onto lines of at most `width` characters, as the paper figures
    write them; past `lines` the rest is cut.
    '''
    wrapped = textwrap.wrap(str(label), width) or [str(label)]
    if len(wrapped) > lines:
        wrapped = wrapped[:lines]
        wrapped[-1] = wrapped[-1][:width - 1].rstrip() + '…'

    return '<br>'.join(wrapped)


def label_width(labels):
    '''Characters of the longest line of the wrapped names.'''
    return max(len(line) for label in labels for line in str(label).split('<br>'))


def row_positions(labels, row):
    '''
    The y of each row, in rows of `row` pixels: a row is one row from the next unless
    their wrapped names need more room, and then the names are as far apart as those of
    two single-line rows, so that a wrapped name reads as one row.
    '''
    line = LABEL_LINE / row
    halves = [line * (str(label).count('<br>') + 1) / 2 for label in labels]
    ys = [0.0]
    for above, below in zip(halves, halves[1:]):
        ys.append(ys[-1] + max(1, above + below + 1 - line))

    return ys


TITLE = 'OrthoHPI 2.0'
SUBTITLE = 'Orthology Prediction of Host-Parasite PPIs'


def navigation():
    '''
    The menu across the top, holding every page; the first one is opened when the
    address names none.
    '''
    pages = [st.Page(path, title=label, icon=icon, default=i == 0)
             for i, (label, icon, path) in enumerate(PAGES)]

    return st.navigation(pages, position='top')


def show_header():
    '''The name of the app and what it does, which the home page opens with.'''
    st.markdown(f"<h1 style='text-align: center; color: #023858;'>{TITLE}</h1>",
                unsafe_allow_html=True)
    st.markdown(f"<h3 style='text-align: center; color: #2b8cbe;'>{SUBTITLE}</h3>",
                unsafe_allow_html=True)

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


def load_config(config_file):
    '''The parsed config, read again only when the file changes.'''
    return read_config(config_file, os.path.getmtime(config_file))


@st.cache_data(show_spinner=False)
def read_config(config_file, modified_at):
    return utils.read_config(config_file)


@st.cache_data(show_spinner=False)
def load_predictions(data_dir):
    '''Every predicted interaction, as pipeline/main.py writes them.'''
    predictions = utils.read_parquet_file(input_file=f'{data_dir}/predictions.parquet')
    predictions['weight'] = predictions['weight'].astype(float)

    return predictions


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


def go_background(data_dir, taxids, niche=None):
    '''
    The GO annotations an enrichment of these species is tested against: those of the
    proteins the network could have been drawn from (filtered_pool).
    '''
    species = [int(t) for t in taxids]
    # fastparquet prunes row groups only, so the exact selection is still applied afterwards
    go_df = utils.read_parquet_file(input_file=f'{data_dir}/gos.parquet',
                                    filters=[('taxid', 'in', species)])
    go_df = go_df[go_df['taxid'].isin(species)]
    pool = filtered_pool(data_dir, tuple(str(t) for t in taxids), niche=niche)
    # a directory carrying neither table is left on the proteome rather than emptied
    if pool:
        go_df = go_df[go_df['#string_protein_id'].isin(pool)]

    return go_df


@st.cache_data(show_spinner=False)
def load_host_tissues(data_dir):
    '''
    Every host protein the TISSUES filter passes with its tissues, before the DeepLoc
    filter, as scripts/build_host_tissues.py writes it; None when it has not been built.
    '''
    tissue_file = os.path.join(data_dir, 'host_tissues.parquet')
    if not os.path.exists(tissue_file):
        return None

    return utils.read_parquet_file(input_file=tissue_file)


def expressed_infected_tissue_proteins(data_dir, config, parasite_taxid):
    '''
    The host proteins annotated to a tissue this parasite infects, whether or not DeepLoc
    later removes them; the filtered pool's when host_tissues.parquet has not been built.
    '''
    tissues = load_host_tissues(data_dir)
    if tissues is None:
        return infected_tissue_proteins(data_dir, config, parasite_taxid)

    mapped = config['tissues']
    infected = {mapped[t].lower() for t in config['parasites'][int(parasite_taxid)]['tissues']}

    return set(tissues.loc[tissues['tissue'].isin(infected), 'protein'])


@st.cache_data(show_spinner=False)
def infected_tissue_proteins(data_dir, config, parasite_taxid):
    '''
    The host proteins annotated to a tissue this parasite is known to infect, over every
    host at once.
    '''
    tissue_file = os.path.join(data_dir, 'tissues_cell_types.parquet')
    if not os.path.exists(tissue_file):
        return set()
    # two of the table's columns, rather than load_tissue_annotation's copy of all of it
    tissues = pd.read_parquet(tissue_file, columns=['Gene', 'Tissue'])

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
def load_protein_sequences(data_dir):
    '''
    Protein sequences, keyed by STRING id, written by pipeline/build_protein_sequences.py.
    '''
    input_file = os.path.join(data_dir, 'protein_sequences.parquet')
    if not os.path.exists(input_file):
        return {}
    sequences = utils.read_parquet_file(input_file=input_file)

    return dict(sequences[['protein', 'sequence']].values)


def to_fasta(records, width=60):
    '''FASTA text of (header, sequence) records, wrapped at width residues.'''
    lines = []
    for header, sequence in records:
        lines.append(f'>{header}')
        lines.extend(sequence[i:i + width] for i in range(0, len(sequence), width))

    return '\n'.join(lines) + '\n'


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
NICHE_ORDER = ['Intracellular', 'Extracellular']
# greys, so the niche is not read as a third categorical split beside the clades and the
# DeepLoc classes
NICHE_COLORS = {'Extracellular': '#c7c7c7', 'Intracellular': '#3d3d3d',
                UNKNOWN_NICHE: '#f0f0f0'}
NICHE_TITLE = 'intracellular / extracellular'


def niche_rank(niche):
    '''Where a niche sorts among the parasites of a group: NICHE_ORDER, then unknown.'''
    return NICHE_ORDER.index(niche) if niche in NICHE_ORDER else len(NICHE_ORDER)


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


# the confidence score is shared across pages the same way as the host
MIN_SCORE, MAX_SCORE, DEFAULT_SCORE = 0.35, 0.9, 0.35
SCORE_STATE_KEY = 'selected_score'
SCORE_WIDGET_KEY = 'selected_score_widget'


def score_slider(help=None):
    '''Draws the confidence score slider, shared across pages through SCORE_STATE_KEY.'''
    if SCORE_WIDGET_KEY not in st.session_state:
        st.session_state[SCORE_WIDGET_KEY] = st.session_state.get(SCORE_STATE_KEY, DEFAULT_SCORE)

    score = st.slider('Confidence score', MIN_SCORE, MAX_SCORE, key=SCORE_WIDGET_KEY, help=help)
    st.session_state[SCORE_STATE_KEY] = score

    return score


def filter_tissues(config, df):
    source = df['taxid1'].unique()[0]
    mapped_tissues = config['tissues']
    tissues = [mapped_tissues[t].lower() for t in config['parasites'][int(source)]['tissues']]
    df = df[df['Tissue'].isin(tissues)]
    
    return df


# minimum expression for a host protein to count in a cell type: HPA's own detection
# threshold, which the mouse and pig atlases share since they are pooled the same way
CELL_TYPE_CUTOFF = 1
# the slider's top; about the 90th percentile of the annotated expression values
CELL_TYPE_MAX_CUTOFF = 200


def keep_expressed_cell_types(df, cutoff=CELL_TYPE_CUTOFF):
    '''The rows of a tissue-annotated frame where a host protein counts as expressed.'''
    return df[df['nCPM'] > cutoff]


def expression_slider(key, help=None):
    '''The nCPM a host protein must exceed in a cell type to count as expressed there.'''
    return st.slider('Cell-type expression (nCPM)', CELL_TYPE_CUTOFF, CELL_TYPE_MAX_CUTOFF,
                     CELL_TYPE_CUTOFF, key=key,
                     help=help or 'A host protein counts towards a cell type when its '
                                  'expression there is above this. 1 nCPM is HPA\'s own '
                                  'detection threshold; the mouse and pig atlases are put on '
                                  'the same scale.')


def expression_phrase(cutoff):
    '''How the captions say what counts as expressed.'''
    return f'above {cutoff:g} nCPM'


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