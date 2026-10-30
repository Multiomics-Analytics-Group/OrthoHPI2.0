import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import utils
import web_utils
import similarity_heatmap
import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()


# fallback for a parasite without a `group` in the config
UNKNOWN_GROUP = 'Unclassified'
UNKNOWN_COLOR = '#999999'
MIN_SCORE = web_utils.MIN_SCORE
# height of one legend row below a figure, in pixels
LEGEND_ROW = 22
# room one legend entry takes, in pixels: marker, padding and the longest name
LEGEND_ENTRY = 40
# room left under the last legend of a dot matrix, in pixels
LEGEND_PAD = 14
# ~6 px a character of the longest name in a dot-matrix legend
DOT_LEGEND_CHAR = 6.0
# room a horizontal legend's title takes above its entries
LEGEND_TITLE_ROW = 18
# smallest the parasite names on a matrix axis are written, in points
SMALLEST_LABEL = 9
# pixels per row of the shared-interactor dot plot
DOT_ROW = 19
DOT_CHROME = 240
# session key counting remounts of the shared-interactor matrix, so the same cell can be
# opened twice running
CELL_NONCE_KEY = 'shared_cell_nonce'
# the tissue selector value meaning no tissue filter
ALL_TISSUES = 'All tissues'


# localization of a host protein, drawn as a band beside its row; NO_LOCALISATION is a
# protein DeepLoc was never run on
NO_LOCALISATION = 'Not available'
LOCALISATION_ORDER = list(web_utils.HOST_CLASSES) + [web_utils.SEVERAL,
                                                     web_utils.NOT_SURFACE, NO_LOCALISATION]
# paler than the grey of a rejected protein
LOCALISATION_COLORS = {**web_utils.LOCALISATION_COLORS, NO_LOCALISATION: '#f0f0f0'}
# names the DeepLoc columns are read under in the hover
DEEPLOC_LABELS = {'surface': 'DeepLoc', 'localizations': 'localizations',
                  **{column: f'P({name.lower()})'
                     for name, column in web_utils.DEEPLOC_SCORES.items()}}


@st.cache_resource(show_spinner=False)
def load_tissue_expression(data_dir):
    '''
    The tissue and cell-type table keyed by host protein, shared rather than copied on
    every call; only ever read.
    '''
    tissues = utils.read_parquet_file(input_file=f'{data_dir}/tissues_cell_types.parquet')

    return tissues.rename({'Gene': 'target'}, axis=1)[['target', 'Tissue', 'Cell type',
                                                       'nCPM']]


def infected_tissues(config):
    '''The tissues each parasite is known to infect (config['parasites']).'''
    mapped_tissues = config['tissues']

    return pd.DataFrame([(str(taxid), mapped_tissues[t].lower())
                         for taxid, parasite in config['parasites'].items()
                         for t in parasite['tissues']],
                        columns=['taxid1', 'Tissue'])


@st.cache_resource(show_spinner=False, max_entries=4)
def get_infected_tissue_expression(data_dir, config, host_taxids):
    '''
    Every prediction against the host, with the cell types its host protein is expressed
    in within the tissues its parasite infects. Joined once a host, at every score, so a
    new threshold only filters it; shared rather than copied, and only ever read.
    '''
    predictions = web_utils.get_host_predictions(data_dir, host_taxids)[
        ['taxid1', 'taxid1_label', 'source', 'target', 'target_name', 'weight']]
    aux = predictions.drop_duplicates().astype({'taxid1': str})
    aux = pd.merge(aux, load_tissue_expression(data_dir), on='target')

    return pd.merge(aux, infected_tissues(config), on=['taxid1', 'Tissue'])


@st.cache_data(show_spinner=False)
def count_interactions_per_tissue(data_dir, config, host_taxids, score=MIN_SCORE,
                                  expression=web_utils.CELL_TYPE_CUTOFF):
    '''
    Predicted interactions per parasite and tissue, and per parasite, tissue and cell
    type, keeping only the tissues each parasite is known to infect
    (config['parasites']).
    '''
    aux = get_infected_tissue_expression(data_dir, config, host_taxids)
    aux = aux[aux['weight'] >= score]

    # keyed by display gene name so aliases of the same host gene do not inflate the counts
    def pairs(frame, *level):
        return (frame.drop_duplicates(list(level) + ['source', 'target_name'])
                     .groupby(list(level), observed=True).size()
                     .rename('interactions').reset_index())

    return (pairs(aux, 'taxid1_label', 'Tissue'),
            pairs(web_utils.keep_expressed_cell_types(aux, expression), 'taxid1_label', 'Tissue',
                  'Cell type'))


def plot_room(labels, column):
    '''
    Pixels of a `column`-wide figure left for the plot itself, once rows named by
    `labels` have taken the room their names need.
    '''
    # 6.2 px a character of the default axis font, 30 px free at either end of the axis
    return column - (6.2 * web_utils.label_width(labels) + 30)


def dot_size(parasites, labels, column):
    '''
    Diameter of the largest dot of a matrix of `parasites` columns whose rows are named
    by `labels`, drawn in a column `column` pixels wide.
    '''
    return min(15, max(5, plot_room(labels, column) / len(parasites)))


def slant_overhang(names, font, width):
    '''
    Pixels the names of the last columns, slanted at 60 degrees, reach past the right
    edge of a plot `width` pixels wide.
    '''
    column = width / len(names)
    # ~0.6 of the font size a character, half of it across at 60 degrees, and the height
    # of the letters
    reach = [0.3 * font * len(name) + font - column * (len(names) - 0.5 - i)
             for i, name in enumerate(names)]

    return max(0, round(max(reach)))


@st.cache_data(show_spinner=False)
def generate_tissue_dots(per_tissue, groups, group_order, niches, palette, column):
    '''
    A dot wherever a parasite is predicted to interact with a host protein expressed in
    a tissue it infects, sized by how many such interactions there are.
    '''
    dots = per_tissue.copy()
    dots['group'] = dots['taxid1_label'].map(lambda p: groups.get(p, UNKNOWN_GROUP))
    dots['parasite'] = dots['taxid1_label'].map(lambda p: f'{p[0]}. {p.split(" ")[1]}')
    dots['niche'] = dots['taxid1_label'].map(niches).fillna(web_utils.UNKNOWN_NICHE)
    order = parasite_order(dots['taxid1_label'].unique(), groups, group_order, niches,
                           ORDER_BY_GROUP)
    parasites = [f'{p[0]}. {p.split(" ")[1]}' for p in order]
    reach = dots.groupby('Tissue').agg(parasites=('taxid1_label', 'nunique'),
                                       total=('interactions', 'sum'))
    tissues = list(reach.sort_values(['parasites', 'total'], ascending=False, kind='stable').index)

    # the axes count cells and the names are ticks, so the strips can be drawn without
    # becoming a category
    dots = dots.assign(x=dots['parasite'].map({p: i for i, p in enumerate(parasites)}),
                       y=dots['Tissue'].map({t: i for i, t in enumerate(tissues)}))

    size = dot_size(parasites, tissues, column)
    sizeref = max(1, dots['interactions'].max()) / size ** 2
    figure = go.Figure()
    # one trace per group, so the groups are the legend
    for group in [g for g in list(palette) + [UNKNOWN_GROUP] if g in set(dots['group'])]:
        rows = dots[dots['group'] == group]
        figure.add_trace(go.Scatter(
            x=rows['x'], y=rows['y'], mode='markers', name=group,
            marker=dict(color=palette.get(group, UNKNOWN_COLOR), size=rows['interactions'],
                        sizemode='area', sizeref=sizeref, sizemin=4, line=dict(width=0)),
            # a dot wider than its column is not cut at either end of the axis
            cliponaxis=False,
            customdata=rows[['Tissue', 'parasite', 'interactions']].to_numpy(),
            hovertemplate='%{customdata[0]}<br>%{customdata[1]}<br>predicted interactions: '
                          '%{customdata[2]}<extra></extra>'))

    groups_shown, niches_shown = add_parasite_strips(figure, dots, parasites, palette)
    font = max(SMALLEST_LABEL, min(11, round(size / 1.1)))
    # the slanted names of the last columns run past the plot, so the right margin holds them
    right = slant_overhang(parasites, font, plot_room(tissues, column))
    height = max(420, 19 * len(tissues) + 240)
    room = add_dot_legend(figure, groups_shown, [], niches_shown, palette,
                          plot_room(tissues, column) - right, height)

    figure.update_layout(height=height, plot_bgcolor='white',
                         # the legends are placed against the foot of the figure, under the
                         # plot; no pad, so the names sit against the strips
                         margin=dict(l=0, r=right, t=10, b=room, pad=0),
                         xaxis_title=None, yaxis_title='tissue the parasite infects')
    # every parasite is named (plotly thins labels that no longer fit); the range runs back
    # past the two strips
    figure.update_xaxes(range=[-0.5, len(parasites) - 0.5], side='top', tickmode='array',
                        tickvals=list(range(len(parasites))), ticktext=parasites,
                        tickangle=-60, automargin=True, ticks='',
                        tickfont=dict(size=font),
                        showgrid=True, gridcolor='#f0f0f0', zeroline=False)
    # reversed, so the tissue the most parasites reach is the top row
    # held off the axis by half the largest dot, which may reach past it
    figure.update_yaxes(range=[len(tissues) - 0.5, -3], tickmode='array',
                        tickvals=list(range(len(tissues))), ticktext=tissues, ticks='',
                        ticklabelstandoff=round(size / 2),
                        automargin=True, showgrid=True, gridcolor='#f0f0f0', zeroline=False)

    return figure


@st.cache_data(show_spinner=False)
def generate_cell_type_bars(per_cell_type, tissue, groups, palette):
    '''
    The cell types of one tissue, each a bar of the interactions predicted there, split
    by taxonomic group.
    '''
    data = per_cell_type[per_cell_type['Tissue'] == tissue].copy()
    data['group'] = data['taxid1_label'].map(lambda p: groups.get(p, UNKNOWN_GROUP))
    data['parasite'] = data['taxid1_label'].map(lambda p: f'{p[0]}. {p.split(" ")[1]}')
    totals = data.groupby('Cell type')['interactions'].sum().sort_values(ascending=False,
                                                                        kind='stable')
    cell_types = list(totals.index)

    figure = px.bar(data, x='interactions', y='Cell type', color='group', orientation='h',
                    color_discrete_map=palette, category_orders={
                        'Cell type': cell_types,
                        'group': [g for g in palette if g in set(data['group'])]},
                    custom_data=['parasite'])
    # two parasites of one group would read as a single bar without a line to part them
    figure.update_traces(marker_line=dict(color='white', width=1),
                         hovertemplate='%{y}<br>%{customdata[0]}<br>predicted interactions: '
                                       '%{x}<extra></extra>')
    figure.update_layout(height=max(320, 26 * len(cell_types) + 160), plot_bgcolor='white',
                         margin=dict(l=0, r=0, t=10, b=10), legend_title_text='',
                         legend=dict(orientation='h', yanchor='bottom', y=1.01, x=0),
                         xaxis_title='predicted interactions', yaxis_title=None,
                         bargap=0.25)
    # the bars are stacked, so the axis has to reach the total
    web_utils.count_ticks(figure, totals.max(), showgrid=True, gridcolor='#f0f0f0')

    return figure

@st.cache_data(show_spinner=False)
def get_tissue_expressed_predictions(data_dir, config, host_taxids, score=MIN_SCORE,
                                     tissue=None):
    '''
    Predictions restricted to the host proteins that are expressed in a tissue the
    parasite is known to infect (config['parasites'][taxid]['tissues']), which is the
    same restriction the tissue matrix applies.
    '''
    predictions = web_utils.get_host_predictions(data_dir, host_taxids)[['taxid1', 'taxid1_label', 'source',
                                                               'target', 'target_name', 'weight']]
    predictions = predictions[predictions['weight'] >= score].drop('weight', axis=1)
    expressed = load_tissue_expression(data_dir)[['target', 'Tissue']].drop_duplicates()

    aux = predictions.drop_duplicates().astype({'taxid1': str})
    aux = pd.merge(aux, expressed, on='target')
    aux = pd.merge(aux, infected_tissues(config), on=['taxid1', 'Tissue'])
    if tissue is not None:
        aux = aux[aux['Tissue'] == tissue]

    return aux[['taxid1_label', 'source', 'target', 'target_name']].drop_duplicates()


# what the matrix axes and dot plot columns can be ordered on; the order decides which
# blocks are visible
ORDER_BY_GROUP = 'taxonomic group'
ORDER_BY_NICHE = 'intracellular / extracellular'


def parasite_order(labels, groups, group_order, niches, order_by):
    '''
    The parasites in the order the figures put them on their axes: the chosen annotation
    first, so its values are contiguous, then the other one, then the name.
    '''
    def rank(parasite):
        clade = group_order.get(groups.get(parasite), len(group_order))
        niche = web_utils.NICHE_ORDER.index(niches[parasite]) \
            if niches.get(parasite) in web_utils.NICHE_ORDER else len(web_utils.NICHE_ORDER)

        return (niche, clade, parasite) if order_by == ORDER_BY_NICHE else (clade, niche,
                                                                           parasite)

    return sorted(labels, key=rank)


@st.cache_data(show_spinner=False)
def get_shared_interactor_similarity(df_pred, groups, group_order, niches, order_by):
    '''
    How alike the host interactors of each pair of parasites are, as the Jaccard
    similarity of the two sets: the host proteins both parasites reach, over the host
    proteins either of them reaches.
    '''
    targets = {g: set(df['target']) for g, df in df_pred.groupby('taxid1_label')}
    targets = {g: t for g, t in targets.items() if t}
    labels = parasite_order(targets, groups, group_order, niches, order_by)
    if len(labels) < 3:
        return None

    shared = np.array([[len(targets[a] & targets[b]) for b in labels] for a in labels],
                      dtype=float)
    # the union of two non-empty sets is non-empty
    union = np.array([[len(targets[a] | targets[b]) for b in labels] for a in labels],
                     dtype=float)
    similarity = shared / union
    np.fill_diagonal(similarity, np.nan)
    np.fill_diagonal(shared, np.nan)

    return (pd.DataFrame(similarity, index=labels, columns=labels),
            pd.DataFrame(shared, index=labels, columns=labels),
            [groups.get(g, UNKNOWN_GROUP) for g in labels],
            [niches.get(g, web_utils.UNKNOWN_NICHE) for g in labels])


# characters a line of a host protein name on the dot matrix is wrapped at
LABEL_WRAP = 30


def label_proteins(df_pred, annotations):
    '''
    Names each host protein by its gene symbol and its descriptive protein name, since
    the symbol on its own identifies the protein only for someone who already knows it.
    '''
    described = {}
    for protein, name in df_pred[['target', 'target_name']].drop_duplicates().values:
        description = (annotations or {}).get(protein, '')
        if description and not description.lower().startswith('uncharacterized'):
            described.setdefault(name, description.rstrip('.'))

    labels = {}
    for name in df_pred['target_name'].unique():
        description = described.get(name)
        labels[name] = f'{name} · {description}' if description else str(name)

    return labels


def summarise_localisations(df_pred, localisations):
    '''
    The DeepLoc call of each host protein of the matrix, keyed by the gene symbol the
    rows are drawn under rather than by the STRING id DeepLoc wrote it against.
    '''
    if localisations is None or localisations.empty:
        return None

    proteins = df_pred[['target', 'target_name']].drop_duplicates()
    called = proteins.merge(localisations, left_on='target', right_on='protein', how='inner')
    if called.empty:
        return None

    scores = [c for c in web_utils.DEEPLOC_SCORES.values() if c in called]
    called['best'] = called[scores].max(axis=1)
    called = called.sort_values('best', ascending=False, kind='stable')
    called = called.drop_duplicates('target_name').set_index('target_name')

    called['surface'] = web_utils.classify_localisation(called)

    return called[['surface'] + scores + ['localizations']]


@st.dialog('Host interactors shared by a pair of parasites', width='large')
def show_shared_interactors_dialog(first, second, df_pred, annotations, localisations):
    '''
    The host proteins behind one cell of the shared-interactor matrix: the cell gives a
    count, and this is what the count is made of.
    '''
    edges = df_pred[df_pred['taxid1_label'].isin([first, second])]
    targets = {p: set(df['target']) for p, df in edges.groupby('taxid1_label')}
    shared = targets.get(first, set()) & targets.get(second, set())

    st.markdown(f'**{first}** and **{second}**', unsafe_allow_html=True)
    st.caption(f'The {len(shared)} host proteins both parasites are predicted to interact '
               f'with, at the confidence the page is set to and among the proteins '
               f'expressed in a tissue each parasite infects. A column counts the proteins '
               f'of that parasite predicted to reach the host protein.')
    if not shared:
        return

    edges = edges[edges['target'].isin(shared)]
    # one column per parasite, indexed by the shared host protein
    degree = edges.groupby(['taxid1_label', 'target'])['source'].nunique().unstack('taxid1_label')
    names = edges.drop_duplicates('target').set_index('target')['target_name']
    labels = label_proteins(edges, annotations)

    table = pd.DataFrame({'Host protein': [labels.get(names[t], names[t]) for t in degree.index],
                          'Identifier': list(degree.index)})
    columns = []
    for parasite in (first, second):
        # named as the matrix axes name a parasite
        column_name = f'{parasite[0]}. {parasite.split(" ")[1]} proteins'
        table[column_name] = degree[parasite].fillna(0).astype(int).values
        columns.append(column_name)

    surface = summarise_localisations(edges, localisations)
    if surface is not None:
        table['DeepLoc'] = [surface['surface'].get(names[t], NO_LOCALISATION)
                            for t in degree.index]

    # the proteins reached by the most parasite proteins first
    table = table.assign(_reach=table[columns].sum(axis=1)).sort_values(
        ['_reach', 'Host protein'], ascending=[False, True], kind='stable').drop(columns='_reach')

    st.dataframe(table, width='stretch', height=400, hide_index=True)


@st.cache_data(show_spinner=False)
def get_top_shared_proteins(df_pred, groups, group_order, niches, order_by,
                            annotations=None, localisations=None, top=40):
    '''The host proteins that the most parasites are predicted to interact with.'''
    edges = df_pred[['taxid1_label', 'source', 'target', 'target_name']].drop_duplicates()
    # one row per dot, keyed by display gene name so aliases do not stack
    pairs = edges[['taxid1_label', 'target_name']].drop_duplicates()
    counts = pairs.groupby('target_name')['taxid1_label'].nunique()
    counts = counts[counts > 1].sort_values(ascending=False, kind='stable')
    if counts.empty:
        return None

    proteins = list(counts.head(top).index)
    degree = edges.groupby(['taxid1_label', 'target_name'])['source'].nunique()
    dots = pairs[pairs['target_name'].isin(proteins)].copy()
    dots['group'] = dots['taxid1_label'].map(lambda p: groups.get(p, UNKNOWN_GROUP))
    dots['parasites'] = dots['target_name'].map(counts)
    dots['degree'] = pd.MultiIndex.from_frame(dots[['taxid1_label', 'target_name']]).map(degree)
    dots['parasite'] = dots['taxid1_label'].map(lambda p: f'{p[0]}. {p.split(" ")[1]}')
    full_labels = label_proteins(df_pred, annotations)
    labels = {name: web_utils.wrap_label(label, LABEL_WRAP)
              for name, label in full_labels.items()}
    dots['protein'] = dots['target_name'].map(labels)
    dots['protein_full'] = dots['target_name'].map(full_labels)
    surface = summarise_localisations(df_pred, localisations)
    if surface is not None:
        dots = dots.join(surface, on='target_name')
        dots['surface'] = dots['surface'].fillna(NO_LOCALISATION)
        dots['localizations'] = dots['localizations'].fillna('')
    order = parasite_order(dots['taxid1_label'].unique(), groups, group_order, niches,
                           order_by)
    dots['niche'] = dots['taxid1_label'].map(niches).fillna(web_utils.UNKNOWN_NICHE)

    return (dots, [labels[p] for p in proteins],
            [f'{p[0]}. {p.split(" ")[1]}' for p in order], len(counts))


@st.cache_data(show_spinner=False)
def get_shared_protein_table(df_pred, annotations=None, localisations=None):
    '''
    Every host protein more than one parasite is predicted to interact with, for
    download, since the dot plot shows only the top of it. One row per protein.
    '''
    edges = df_pred[['taxid1_label', 'source', 'target', 'target_name']].drop_duplicates()
    degree = edges.groupby(['target_name', 'taxid1_label'])['source'].nunique()
    counts = degree.groupby(level='target_name').size()
    counts = counts[counts > 1].sort_values(ascending=False, kind='stable')
    ids = edges.groupby('target_name')['target'].agg(lambda t: ', '.join(sorted(set(t))))
    labels = label_proteins(df_pred, annotations)

    table = pd.DataFrame({
        'Host protein': counts.index,
        'STRING id': [ids[p] for p in counts.index],
        'Description': [labels[p].partition(' · ')[2] for p in counts.index],
        'Parasites': counts.values,
        'Parasites (interacting parasite proteins)': [
            '; '.join(f'{parasite} ({n})' for parasite, n in degree[p].items())
            for p in counts.index]})
    surface = summarise_localisations(df_pred, localisations)
    if surface is not None:
        table['DeepLoc'] = [surface['surface'].get(p, NO_LOCALISATION) for p in counts.index]

    return table


def localisation_labels(dots, proteins):
    '''
    What the band beside a row of the dot plot says when it is hovered: the class of the
    host protein, and for a protein called for more than one, which classes those are.
    '''
    called = dots.drop_duplicates('protein').set_index('protein')
    scores = {c: web_utils.DEEPLOC_SCORES[c] for c in web_utils.HOST_CLASSES
              if web_utils.DEEPLOC_SCORES[c] in called}
    labels = []
    for protein in proteins:
        row = called.loc[protein]
        crossed = [c for c, column in scores.items()
                   if row[column] > web_utils.DEEPLOC_CUTOFFS[c]]
        # a protein of one class is named by it, one with no localisation by what the legend
        # calls it
        labels.append(' + '.join(crossed) if len(crossed) > 1 else row['surface'])

    return labels


def add_parasite_strips(figure, dots, parasites, palette):
    '''The taxonomic group and the niche of each parasite, as two bands above its column.'''
    clade_of = dict(zip(dots['parasite'], dots['group']))
    niche_of = dict(zip(dots['parasite'], dots['niche']))
    clades = [clade_of.get(p, UNKNOWN_GROUP) for p in parasites]
    niches = [niche_of.get(p, web_utils.UNKNOWN_NICHE) for p in parasites]
    groups_shown = [g for g in list(palette) + [UNKNOWN_GROUP] if g in set(clades)]
    niches_shown = [n for n in web_utils.NICHE_ORDER + [web_utils.UNKNOWN_NICHE]
                    if n in set(niches)]

    # the group a cell clear of the first row and the niche outside it
    for offset, (values, colors) in enumerate(
            [(clades, {g: palette.get(g, UNKNOWN_COLOR) for g in groups_shown}),
             (niches, {n: web_utils.NICHE_COLORS[n] for n in niches_shown})]):
        scale, codes = similarity_heatmap.flat_colour_scale(values, colors)
        figure.add_trace(go.Heatmap(z=[codes], x=list(range(len(parasites))),
                                    y0=-1.4 - offset, dy=1, text=[values], xgap=1, **scale))

    return groups_shown, niches_shown


def strip_width(parasites, rows, height, room):
    '''
    Width in x units of the band beside the rows of the dot plot, at which its cells
    come out as wide as they are high; `rows` is the span of the rows in y units.
    '''
    # the y axis is opened up past the top row to hold the two strips
    row = (height - DOT_CHROME) / (rows + 2.5)
    # room = row * (parasites + 1.5 * width) / width, solved for the width
    return row * parasites / max(room - 1.5 * row, row)


def add_localisation_strip(figure, dots, proteins, ys, width):
    '''Where DeepLoc puts each host protein, as a band beside its row.'''
    surface_of = dict(zip(dots['protein'], dots['surface']))
    localisations = [surface_of.get(p, NO_LOCALISATION) for p in proteins]
    shown = [c for c in LOCALISATION_ORDER if c in set(localisations)]
    scale, codes = similarity_heatmap.flat_colour_scale(localisations,
                                     {c: LOCALISATION_COLORS[c] for c in shown})
    # the band is `width` x units across, an x unit being a parasite
    figure.add_trace(go.Heatmap(z=[[code] for code in codes], x0=-0.5 - width, dx=width,
                                y=ys,
                                text=[[label] for label in
                                      localisation_labels(dots, proteins)],
                                ygap=1, **scale))

    return shown


def add_dot_legend(figure, groups, localisations, niches, palette, width, height):
    '''
    Names the two strips the dot plot carries no legend for, and lays the three keys of
    the figure out under it as three legends of their own.
    '''
    # dots are sized by degree, so the legend is drawn from markers of one size instead of
    # the traces
    for trace in figure.data:
        trace.update(legendgroup=trace.name, showlegend=False)
    for name, color in [(g, palette.get(g, UNKNOWN_COLOR)) for g in groups]:
        figure.add_scatter(x=[None], y=[None], mode='markers', name=name, legendgroup=name,
                           marker=dict(symbol='circle', size=10, color=color),
                           hoverinfo='skip', showlegend=True)
    for name, color, legend in ([(c, LOCALISATION_COLORS[c], 'legend2') for c in localisations]
                                + [(n, web_utils.NICHE_COLORS[n], 'legend3') for n in niches]):
        figure.add_scatter(x=[None], y=[None], mode='markers', name=name, legend=legend,
                           marker=dict(symbol='square', size=10, color=color),
                           hoverinfo='skip', showlegend=True)

    def pixels(names):
        '''The room a legend takes: its title, and the rows its entries wrap onto.'''
        entry = LEGEND_ENTRY + DOT_LEGEND_CHAR * max(len(name) for name in names)

        return LEGEND_TITLE_ROW + LEGEND_ROW * -(-len(names)
                                                 // max(1, int(width // entry)))

    keys = [(legend, title, names, pixels(names))
            for legend, title, names in [('legend', 'taxonomic group', groups),
                                         ('legend2', 'DeepLoc', localisations),
                                         ('legend3', web_utils.NICHE_TITLE, niches)]
            if names]
    room = sum(taken for *_, taken in keys) + LEGEND_PAD

    # measured from the foot of the figure, since the top margin plotly takes for the names
    # is not known here
    offset = room / height
    for legend, title, names, taken in keys:
        figure.update_layout({legend: dict(
            orientation='h', yanchor='top', y=offset, yref='container', x=0,
            traceorder='normal', title_text=title, title_font=dict(size=11),
            font=dict(size=11))})
        offset -= taken / height

    return room


@st.cache_data(show_spinner=False)
def generate_shared_protein_dots(dots, proteins, parasites, palette, column):
    '''
    A dot wherever a parasite is predicted to interact with one of the proteins, the
    parasites in the order the other figures use so the taxonomic groups stay together,
    and the proteins ordered by how many parasites reach them.
    '''
    localised = 'surface' in dots.columns
    # a wrapped name opens its row up
    ys = web_utils.row_positions(proteins, DOT_ROW)
    cells = {'x': {p: i for i, p in enumerate(parasites)},
             'y': dict(zip(proteins, ys))}
    dots = dots.assign(x=dots['parasite'].map(cells['x']),
                       y=dots['protein'].map(cells['y'])).dropna(subset=['x', 'y'])

    # written out: `parasites` and `degree` are two different counts, and the parasite is
    # named from the row rather than the x
    hover_columns = ['protein_full', 'parasites', 'degree', 'parasite']
    hover_lines = ['%{customdata[0]}', 'parasites reaching it: %{customdata[1]}',
                   'proteins of %{customdata[3]} reaching it: %{customdata[2]}']
    if localised:
        # the probabilities behind the band, two to a line
        scores = [web_utils.DEEPLOC_SCORES[c] for c in web_utils.HOST_CLASSES
                  if web_utils.DEEPLOC_SCORES[c] in dots]
        parts = [f'{DEEPLOC_LABELS[name]} %{{customdata[{len(hover_columns) + i}]:.2f}}'
                 for i, name in enumerate(scores)]
        hover_columns += scores
        hover_lines += [', '.join(parts[i:i + 2]) for i in range(0, len(parts), 2)]

    size = dot_size(parasites, proteins, column)
    # the area of the largest dot stands for the largest degree
    sizeref = max(1, dots['degree'].max()) / size ** 2
    figure = go.Figure()
    # one trace per group, so the groups are the legend
    for group in [g for g in list(palette) + [UNKNOWN_GROUP] if g in set(dots['group'])]:
        rows = dots[dots['group'] == group]
        figure.add_trace(go.Scatter(
            x=rows['x'], y=rows['y'], mode='markers', name=group,
            marker=dict(color=palette.get(group, UNKNOWN_COLOR), size=rows['degree'],
                        sizemode='area', sizeref=sizeref, sizemin=4, line=dict(width=0)),
            customdata=rows[hover_columns].to_numpy(),
            hovertemplate='<br>'.join(hover_lines) + '<extra></extra>'))

    groups_shown, niches = add_parasite_strips(figure, dots, parasites, palette)
    font = max(SMALLEST_LABEL, min(11, round(size / 1.1)))
    # the slanted names of the last columns run past the plot, so the right margin holds them
    right = slant_overhang(parasites, font, plot_room(proteins, column))
    # the legends wrap inside what the row labels leave of the column
    height = max(420, round(DOT_ROW * (ys[-1] + 1)) + DOT_CHROME)
    room = plot_room(proteins, column) - right
    band = strip_width(len(parasites), ys[-1] + 1, height, room) if localised else 0
    localisations = (add_localisation_strip(figure, dots, proteins, ys, band) if localised
                     else [])
    legends = add_dot_legend(figure, groups_shown, localisations, niches, palette, room,
                             height)

    figure.update_layout(height=height, plot_bgcolor='white',
                         # the legends are placed against the foot of the figure, under the
                         # plot; no pad, so the names sit against the strips
                         margin=dict(l=0, r=right, t=10, b=legends, pad=0),
                         xaxis_title=None, yaxis_title='host protein')
    # every parasite is named at the size a column has room for; the range runs back past
    # the band beside the rows
    figure.update_xaxes(range=[-0.5 - 1.5 * band, len(parasites) - 0.5],
                        side='top', tickmode='array', tickvals=list(range(len(parasites))),
                        ticktext=parasites, tickangle=-60, automargin=True, ticks='',
                        tickfont=dict(size=font),
                        showgrid=True, gridcolor='#f0f0f0', zeroline=False)
    # reversed, so the most-shared protein is the top row
    figure.update_yaxes(range=[ys[-1] + 0.5, -3], tickmode='array',
                        tickvals=ys, ticktext=proteins, ticks='',
                        automargin=True, showgrid=True, gridcolor='#f0f0f0', zeroline=False)

    return figure


# the processes the GO dot matrix shows: those enriched for the most parasites, or the
# strongest few of each parasite
GO_ROWS_SHARED = 'enriched for the most parasites'
GO_ROWS_STRONGEST = 'strongest for each parasite'
GO_SHARED_ROWS = 40
GO_STRONGEST_PER_PARASITE = 2
# the colour ramp of the network page's enrichment figures, for the odds ratio here
GO_SEQUENTIAL = ['#bcdcec', '#7fc0dd', '#3f9fca', '#2b8cbe', '#12587d', '#08324a']
# characters of a process name on the axis; the hover has the whole of it
GO_LABEL_WIDTH = 48


@st.cache_resource(show_spinner=False, max_entries=8)
def get_go_background(data_dir, host_taxids, niche):
    '''
    The GO background of a niche, with what the enrichment reads from it worked out; the
    same at every threshold, and only ever read.
    '''
    go_df = web_utils.go_background(data_dir, host_taxids, niche)

    return go_df, utils.prepare_go_background(go_df)


@st.cache_data(show_spinner=False)
def get_parasite_enrichments(df_pred, data_dir, host_taxids, niches):
    '''
    The biological processes over-represented among the host proteins of each parasite,
    each tested against the host proteins of that parasite's niche, as the network page
    tests one.
    '''
    enrichments = []
    for parasite, rows in df_pred.groupby('taxid1_label'):
        go_df, background = get_go_background(data_dir, host_taxids, niches.get(parasite))
        enrichment = utils.calculate_enrichment(set(rows['target']), go_df,
                                                background=background)
        if not enrichment.empty:
            enrichments.append(enrichment.drop(columns=['B', 'C', 'D'])
                               .rename(columns={'A': 'n_proteins'})
                               .assign(parasite=parasite))
    if not enrichments:
        return None

    return pd.concat(enrichments, ignore_index=True)


def select_go_rows(enriched, rows):
    '''
    The processes the dot matrix shows, most widely enriched first: the ones enriched for
    the most parasites, or the union of each parasite's most significant ones.
    '''
    reach = enriched.groupby('go_id').agg(parasites=('parasite', 'nunique'),
                                          best=('fdr_bh', 'min'))
    if rows == GO_ROWS_STRONGEST:
        strongest = (enriched.sort_values(['fdr_bh', 'go_id'])
                     .groupby('parasite').head(GO_STRONGEST_PER_PARASITE)['go_id'])
        reach = reach.loc[strongest.unique()]
    order = reach.sort_values(['parasites', 'best'], ascending=[False, True], kind='stable')

    return list(order.index if rows == GO_ROWS_STRONGEST else order.index[:GO_SHARED_ROWS])


def short_parasite(label):
    return f'{label[0]}. {label.split(" ")[1]}'


@st.cache_data(show_spinner=False)
def generate_go_dots(enriched, terms, parasites, groups, niches, palette, width):
    '''
    A dot wherever a process is enriched among the host proteins of a parasite, sized by
    the number of them annotated to it and shaded by its odds ratio, the parasites in the
    order the other figures use.
    '''
    names = dict(zip(enriched['go_id'], enriched['go_term']))
    labels = [names[t] if len(names[t]) <= GO_LABEL_WIDTH
              else names[t][:GO_LABEL_WIDTH - 1].rstrip() + '…' for t in terms]
    columns = [short_parasite(p) for p in parasites]
    dots = enriched[enriched['go_id'].isin(terms)].copy()
    dots['x'] = dots['parasite'].map({p: i for i, p in enumerate(parasites)})
    dots['y'] = dots['go_id'].map({t: i for i, t in enumerate(terms)})
    odds = pd.to_numeric(dots['odds_ratio'], errors='coerce')
    finite = odds[np.isfinite(odds)]
    # an infinite ratio is drawn at the largest finite one and said so in the hover
    cap = finite.max() if not finite.empty else 1.0
    dots['log_odds'] = np.log2(odds.where(np.isfinite(odds), cap).clip(lower=1))
    dots['odds_text'] = np.where(np.isfinite(odds), odds.map('{:.1f}'.format), 'infinite')
    dots['parasite_name'] = dots['parasite']

    size = dot_size(columns, labels, width)
    sizeref = max(1, dots['n_proteins'].max()) / size ** 2
    figure = go.Figure(go.Scatter(
        x=dots['x'], y=dots['y'], mode='markers', showlegend=False,
        marker=dict(size=dots['n_proteins'], sizemode='area', sizeref=sizeref, sizemin=3,
                    color=dots['log_odds'], colorscale=GO_SEQUENTIAL, line=dict(width=0),
                    colorbar=dict(title='log<sub>2</sub> odds<br>ratio', thickness=12,
                                  outlinewidth=0, len=0.4, y=1, yanchor='top')),
        customdata=dots[['go_term', 'parasite_name', 'n_proteins', 'odds_text',
                         'fdr_bh']].to_numpy(),
        hovertemplate='<b>%{customdata[0]}</b><br>%{customdata[1]}<br>'
                      'host proteins annotated to it: %{customdata[2]}<br>'
                      'odds ratio: %{customdata[3]}<br>FDR: %{customdata[4]:.2e}'
                      '<extra></extra>'))

    strips = pd.DataFrame({'parasite': columns,
                           'group': [groups.get(p, UNKNOWN_GROUP) for p in parasites],
                           'niche': [niches.get(p, web_utils.UNKNOWN_NICHE)
                                     for p in parasites]})
    groups_shown, niches_shown = add_parasite_strips(figure, strips, columns, palette)
    font = max(SMALLEST_LABEL, min(11, round(size / 1.1)))
    # the slanted names of the last columns run past the plot, so the right margin holds them
    right = slant_overhang(columns, font, plot_room(labels, width))
    height = max(420, DOT_ROW * len(terms) + DOT_CHROME)
    legends = add_dot_legend(figure, groups_shown, [], niches_shown, palette,
                             plot_room(labels, width) - right, height)

    figure.update_layout(height=height, plot_bgcolor='white',
                         # no pad, so the names sit against the strips
                         margin=dict(l=0, r=right, t=10, b=legends, pad=0),
                         xaxis_title=None, yaxis_title='biological process')
    figure.update_xaxes(range=[-0.5, len(columns) - 0.5], side='top', tickmode='array',
                        tickvals=list(range(len(columns))), ticktext=columns,
                        tickangle=-60, automargin=True, ticks='',
                        tickfont=dict(size=font),
                        showgrid=True, gridcolor='#f0f0f0', zeroline=False)
    # reversed, so the most widely enriched process is the top row
    figure.update_yaxes(range=[len(terms) - 0.5, -3], tickmode='array',
                        tickvals=list(range(len(terms))), ticktext=labels, ticks='',
                        automargin=True, showgrid=True, gridcolor='#f0f0f0', zeroline=False)

    return figure


# the matrix keeps its cells square, so it needs the column width; drawn for
# DEFAULT_PAGE_WIDTH until the browser answers
column = web_utils.column_width(2)

st.caption('Compare all the parasites of a selected host.')
st.markdown("---")

col1, col2, col3 = st.columns(3)

with col1:
    st.write('')

with col2:
    selected_host, selected_taxids = web_utils.host_selector(
        config, web_utils.load_predictions(data_dir),
        'Select a host')

with col3:
    st.write('')


if selected_host != web_utils.NO_HOST:
    # the heatmap and dot matrix read the same filtered interactions
    parasite_groups = {p['label']: p.get('group', UNKNOWN_GROUP)
                       for p in config['parasites'].values()}
    group_order = {g: i for i, g in enumerate(config.get('parasite_groups', {}))}
    niches = web_utils.get_niches(config)

    # one slider for the four figures below; the tissue dots at the foot keep every
    # prediction
    slider_column, order_column, _ = st.columns(3)
    with slider_column:
        score = web_utils.score_slider(
            help='Interactions predicted below this confidence are left out of '
                 'the four figures below. The tissue plot at the foot of the '
                 'page counts every prediction.')
    with order_column:
        # the order decides which blocks the figures show: contiguous values make a square
        # against the diagonal
        order_by = st.radio('Order the parasites by', [ORDER_BY_NICHE, ORDER_BY_GROUP],
                            horizontal=True,
                            help='Which annotation the two figures below put next to each '
                                 'other on their axes. The strips beside the axes show both '
                                 'either way; this is which of them comes out in blocks.')
    counted = get_tissue_expressed_predictions(data_dir, config, selected_taxids, score)
    shared_similarity = get_shared_interactor_similarity(counted, parasite_groups,
                                                        group_order, niches, order_by)

    per_tissue, per_cell_type = count_interactions_per_tissue(data_dir, config,
                                                              selected_taxids, score)
    ranked = per_tissue.groupby('Tissue')['interactions'].sum().sort_values(ascending=False,
                                                                           kind='stable')
    # the tissues with cell-type data at all, before the expression slider narrows them
    annotated = set(per_cell_type['Tissue'])
    # every tissue is offered, so a missing cell-type annotation is reported rather than
    # guessed at
    choices = list(ranked.index)

    matrix, shared = st.columns(2)

    with matrix:
        st.subheader(
            "Overlap of the host interactors",
            help='How alike the host interactors of each pair of parasites are, as the '
                 'Jaccard similarity of the two sets: the host proteins both reach, over '
                 'the host proteins either of them reaches. A count of shared proteins on '
                 'its own follows how many interactors the pair have between them, and '
                 'ranks the best-predicted parasites above the most alike; this does not. '
                 'Two strips '
                 'run along each axis: the taxonomic group of the parasite, and whether it '
                 'lives inside a host cell or outside one. Whichever the parasites are '
                 'ordered by comes out in blocks. The diagonal, where a parasite meets '
                 'itself, is greyed out. Hover a cell for the count behind its ratio, or '
                 'click it to see the shared host proteins.')
        st.caption('How alike the host proteins targeted by each pair of parasites are. '
                   'Click a cell to see the shared proteins of the pair.')
        if shared_similarity is not None:
            # the figure carries the column width rather than being stretched, which keeps
            # the matrix square after full screen
            similarity, shared_counts, clades, row_niches = shared_similarity
            figure, click_layer = similarity_heatmap.generate(
                similarity, shared_counts,
                [(clades, config.get('parasite_groups', {})),
                 (row_niches, {n: web_utils.NICHE_COLORS[n]
                               for n in web_utils.NICHE_ORDER + [web_utils.UNKNOWN_NICHE]})],
                column)
            # remounted after every dialog: Streamlit drops an identical selection and
            # plotly reads a second click as a deselection
            nonce = st.session_state.get(CELL_NONCE_KEY, 0)
            clicked = st.plotly_chart(
                figure, width='content', on_select='rerun', selection_mode='points',
                key=f'shared_cells_{selected_host}_{order_by}_{score}_{nonce}')
            # the layer the click came from is checked
            cell = next((point for point
                         in (clicked or {}).get('selection', {}).get('points', [])
                         if point.get('curve_number') == click_layer), None)
            if cell is not None:
                # the two parasites are the row and column the marker sits on
                parasites = list(shared_similarity[0].index)
                st.session_state[CELL_NONCE_KEY] = nonce + 1
                show_shared_interactors_dialog(
                    parasites[int(cell['y'])], parasites[int(cell['x'])], counted,
                    web_utils.load_protein_annotations(data_dir),
                    web_utils.load_deeploc_localisations(data_dir))
        else:
            st.text(f'Fewer than three parasites of {selected_host} share any host protein')

    with shared:
        # the description counts what the tissue below leaves but sits above it
        heading = st.empty()
        shared_help = ("A dot wherever a parasite is predicted to interact with one, sized "
                       "by the number of that parasite's proteins reaching it. The band "
                       'beside each row gives the DeepLoc 2 localization of the host '
                       'protein; hover a dot for the probabilities behind it. Above the '
                       'columns run the taxonomic group of the parasite and whether it lives '
                       'inside a host cell or outside one.')
        # offered whatever the tissue leaves, so an emptying choice can be changed back
        shared_tissue = st.selectbox(
            'Tissue', [ALL_TISSUES] + choices, index=0, key='shared_tissue',
            help='Keep only the host proteins expressed in this tissue, for the parasites '
                 'known to infect it. Tissues the parasites infect, most interactions first.')
        # narrowed to one tissue where one is chosen; the heatmap keeps every tissue
        shared_predictions = (counted if shared_tissue == ALL_TISSUES else
                              get_tissue_expressed_predictions(data_dir, config,
                                                               selected_taxids, score,
                                                               shared_tissue))
        top_shared = get_top_shared_proteins(shared_predictions, parasite_groups, group_order,
                                             niches, order_by,
                                             web_utils.load_protein_annotations(data_dir),
                                             web_utils.load_deeploc_localisations(data_dir))
        where = '' if shared_tissue == ALL_TISSUES else f' in {shared_tissue}'
        if top_shared is not None:
            # the count of what the rows were taken from belongs to the caption
            *figure_arguments, shareable = top_shared
            shown = len(figure_arguments[1])
            # a truncated figure says what it is a top of
            if shareable > shown:
                shared_help = (f'The {shown} host proteins{where} reached by the most '
                               f'parasites are shown, of the {shareable} reached by more '
                               'than one. ' + shared_help)
                selection = f'The {shown} host proteins{where} targeted by the most parasites.'
            else:
                selection = (f'The {shown} host proteins{where} targeted by more than one '
                             'parasite, most first.')
            heading.subheader("Most common host interactors", help=shared_help)
            st.caption(selection)
            st.plotly_chart(
                generate_shared_protein_dots(*figure_arguments,
                                             config.get('parasite_groups', {}), column),
                width='stretch')
            shared_table = get_shared_protein_table(
                shared_predictions, web_utils.load_protein_annotations(data_dir),
                web_utils.load_deeploc_localisations(data_dir))
            tissue_suffix = '' if shared_tissue == ALL_TISSUES else f'_{shared_tissue}'
            st.download_button(f'Download all {shareable} shared host interactors',
                               data=utils.convert_df(shared_table),
                               file_name=f'{selected_host}_shared_host_interactors'
                                         f'{tissue_suffix}.tsv',
                               mime='text/csv')
        else:
            heading.subheader("Most common host interactors", help=shared_help)
            st.info(f'No host protein{where} is reached by more than one parasite at this '
                    'confidence.')

    tissues, cell_types = st.columns(2)
    with tissues:
        st.subheader(
            "Interactions per tissue",
            help='Predicted interactions per parasite and tissue, restricted to the tissues '
                 'each parasite is known to infect and sized by the number of interactions '
                 'with proteins expressed there. Tissues are ordered by the number of '
                 'parasites infecting them. Above the columns run the taxonomic group of the '
                 'parasite and whether it lives inside a host cell or outside one, the same '
                 'two strips the figures above carry. An interaction is counted once per '
                 'tissue, irrespective of the number of cell types the host protein is '
                 'expressed in.')
        st.caption('How many interactions each parasite has in each tissue it infects.')
        if per_tissue.empty:
            st.info('No predicted interaction is left at this confidence in a tissue the '
                    'parasites are known to infect.')
        else:
            st.plotly_chart(generate_tissue_dots(per_tissue, parasite_groups, group_order,
                                                 niches, config.get('parasite_groups', {}),
                                                 column),
                            width='stretch')

    with cell_types:
        if choices:
            # the description follows the expression slider but sits above it
            heading = st.empty()
            st.caption('How many interactions take place in each cell type of a tissue.')
            tissue_column, expression_column = st.columns(2)
            with tissue_column:
                tissue = st.selectbox('Tissue', choices, index=0,
                                      help='Tissues the parasites infect, most interactions '
                                           'first')
            with expression_column:
                expression = web_utils.expression_slider('compare_expression')
            # written above the selectors once the slider has been read
            heading.subheader(
                "Interactions per cell type",
                help='Predicted interactions per cell type of the selected tissue, stacked '
                     'by taxonomic group. A host protein counts towards a cell type where '
                     f'it is expressed, {web_utils.expression_phrase(expression)}. A '
                     'protein expressed in several cell types counts in each, so the bars '
                     'are not a partition of the tissue.')
            _, expressed = count_interactions_per_tissue(data_dir, config, selected_taxids,
                                                         score, expression)
            if tissue in annotated and tissue not in set(expressed['Tissue']):
                st.info(f'No host protein interacting with these parasites is expressed '
                        f'{web_utils.expression_phrase(expression)} in a cell type of {tissue}.')
            elif tissue in annotated:
                st.plotly_chart(generate_cell_type_bars(expressed, tissue, parasite_groups,
                                                        config.get('parasite_groups', {})),
                                width='stretch')
            else:
                st.info(f'No single cell data available for {tissue} in {selected_host}, so '
                        'the interactions there cannot be split by cell type.')

    st.subheader(
        "Enriched biological processes",
        help='Gene Ontology biological processes over-represented among the host proteins '
             'each parasite is predicted to reach, tested one parasite at a time as on the '
             "network page: one-sided Fisher's exact test against the host proteins of the "
             "parasite's niche the pipeline passed, corrected with Benjamini-Hochberg. A "
             'dot wherever a process passes the FDR for a parasite, sized by the number of '
             'its host interactors annotated to the process and shaded by the odds ratio. '
             'The odds ratio rather than the FDR is compared, since the FDR follows how '
             'many interactors a parasite has. Intracellular and extracellular parasites '
             'are tested against different host proteins, so a missing dot can come from '
             'the background as well as the biology.')
    st.caption('Which biological processes are over-represented among the host proteins '
               'each parasite targets.')
    rows_column, fdr_column = st.columns([2, 1])
    with rows_column:
        go_rows = st.radio('Processes shown', [GO_ROWS_STRONGEST, GO_ROWS_SHARED],
                           horizontal=True,
                           help=f'The {GO_STRONGEST_PER_PARASITE} most significant processes '
                                'of each parasite, which brings out what sets one apart, or '
                                f'the {GO_SHARED_ROWS} enriched for the most parasites, which '
                                'favours broad processes many share.')
    with fdr_column:
        go_fdr = st.radio('False discovery rate', (0.01, 0.05, 0.1), index=1,
                          horizontal=True, key='compare_go_fdr',
                          help='The Benjamini-Hochberg corrected significance a process has '
                               'to reach to be drawn for a parasite.')
    with st.spinner('Testing the host interactors of every parasite'):
        go_enrichments = get_parasite_enrichments(counted, data_dir, tuple(selected_taxids),
                                                  niches)
    go_enriched = (None if go_enrichments is None else
                   go_enrichments[go_enrichments['fdr_bh'] <= go_fdr])
    if go_enriched is None or go_enriched.empty:
        st.info(f'No biological process passes an FDR of {go_fdr} for any parasite at this '
                'confidence.')
    else:
        go_parasites = parasite_order(counted['taxid1_label'].unique(), parasite_groups,
                                      group_order, niches, order_by)
        st.plotly_chart(generate_go_dots(go_enriched, select_go_rows(go_enriched, go_rows),
                                         go_parasites, parasite_groups, niches,
                                         config.get('parasite_groups', {}),
                                         2 * column + web_utils.COLUMN_GAP),
                        width='stretch')
        st.download_button('Download the enriched processes of every parasite',
                           data=utils.convert_df(go_enriched.drop(columns=['nodes'])),
                           file_name=f'{selected_host}_parasite_go_enrichment.tsv',
                           mime='text/csv')

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()
