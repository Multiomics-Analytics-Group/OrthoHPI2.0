import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import utils
import web_utils
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import body_figure
from css import style

st.set_page_config(layout="wide", page_title="OrthoHPI 2.0", menu_items={})
style.load_css()
web_utils.show_header('Multi-host parasites')

config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()


# opens at the bottom of the range: the multi-host parasites are the small interactomes, and
# 0.7 leaves two of them on the page
MIN_SCORE, MAX_SCORE, DEFAULT_SCORE = 0.35, 0.9, 0.35
# Okabe-Ito without the blues, which say how far a link carried over; the config colours
# give human the grey of an absent host. A host without an entry falls back to its config
# colour
HOST_COLORS = {'Homo sapiens (human)': '#D55E00', 'Sus scrofa (pig)': '#CC79A7',
               'Mus musculus (mouse)': '#009E73', 'Rattus norvegicus (rat)': '#E69F00'}
# a link shared by every host, and the shades a link shared by some of three or more hosts
# is drawn in
SHARED_COLOR = '#08519c'
PARTIAL_COLORS = ['#6baed6', '#9ecae1', '#4292c6', '#c6dbef']
# how far toward white each host after the first sharing a config colour is mixed (rat and
# mouse share one)
SPECIES_TINT = 0.45
# a host the combination does not include; no host is drawn in a grey
ABSENT_COLOR = '#e0e0e0'
# why a host missed a link, in the order explain_host_specific tests them; the paragraph
# under the table counts the rows of each
ABSENT_FAMILY = 'no protein of the family in this host'
NOT_EXPRESSED = 'family present, but no protein of it reaches the infected tissues'
OUT_OF_REACH = 'family expressed, but not where this parasite can reach it'
NOT_TRANSFERRED = 'protein available, link not transferred'
# most gene symbols in an orthology group label before the rest are left to the hover
SYMBOLS_IN_LABEL = 3
# the label is cut here whatever the count, since parasite proteins without a symbol are
# drawn under their locus
LABEL_CHARS = 24
def short_name(parasite):
    '''`Trichinella spiralis` as `T. spiralis`, the abbreviation the other pages use.'''
    parts = parasite.split(' ')

    return f'{parasite[0]}. {parts[1]}' if len(parts) > 1 else parasite


@st.cache_data(show_spinner=False)
def load_host_orthologs(data_dir):
    '''
    Which proteins of each host belong to each of the orthology groups the predictions
    reach, written by scripts/build_host_orthologs.py.
    '''
    input_file = os.path.join(data_dir, 'host_orthologs.parquet')
    if not os.path.exists(input_file):
        return None

    return utils.read_parquet_file(input_file=input_file)


def host_label(host, config):
    '''The configured species name of a host taxid.'''
    return config['hosts'][int(host)]['label']


def host_color(host_label, config):
    '''
    The colour this page fixes for the host (HOST_COLORS), falling back to the one the
    config gives it -- a host added to the config and not to HOST_COLORS still gets
    drawn in the colour the rest of the app knows it by.
    '''
    if host_label in HOST_COLORS:
        return HOST_COLORS[host_label]

    for taxid, host in config['hosts'].items():
        if host['label'] == host_label:
            return host['color']

    return SHARED_COLOR


@st.cache_data(show_spinner=False)
def get_multi_host_predictions(data_dir, config, score=MIN_SCORE):
    '''
    Every prediction of the parasites that are predicted against more than one host,
    labelled with the host it was predicted against.
    '''
    predictions = web_utils.load_predictions(data_dir)
    predictions = predictions[predictions['weight'] >= score].copy()
    predictions['host'] = predictions['taxid2'].map(lambda t: host_label(t, config))

    hosts = predictions.groupby('taxid1')['host'].nunique()

    return predictions[predictions['taxid1'].isin(hosts[hosts > 1].index)]


def combination_label(hosts):
    '''
    The name of a set of hosts: the hosts abbreviated and joined, in a fixed order so
    that the same set is always the same label and the same column of the figures.
    '''
    return ' + '.join(short_name(host) for host in sorted(hosts))


def tint(color, amount):
    '''
    Mixes a colour toward white, the same way the network page lightens a species colour
    into a node fill.
    '''
    color = str(color)
    if not color.startswith('#') or len(color) != 7:
        return color
    channels = [int(color[i:i + 2], 16) for i in (1, 3, 5)]

    return '#%02x%02x%02x' % tuple(round(c + (255 - c) * amount) for c in channels)


def host_palette(all_hosts, config):
    '''One colour per host of the parasite.'''
    sharing = {}
    for host in all_hosts:
        sharing.setdefault(host_color(host, config), []).append(host)

    return {host: tint(color, min(SPECIES_TINT * i, 0.7))
            for color, hosts in sharing.items() for i, host in enumerate(hosts)}


def combination_palette(links, all_hosts, config):
    '''
    The colour each set of hosts is drawn in, over both figures, so that a set is the
    same colour wherever it is read.
    '''
    members = links.drop_duplicates('combination').set_index('combination')['hosts']
    hosts_palette = host_palette(all_hosts, config)
    palette, partial = {}, 0
    for combination in order_combinations(links):
        hosts = members[combination]
        if len(hosts) == len(all_hosts):
            palette[combination] = SHARED_COLOR
        elif len(hosts) == 1:
            palette[combination] = hosts_palette[next(iter(hosts))]
        else:
            palette[combination] = PARTIAL_COLORS[partial % len(PARTIAL_COLORS)]
            partial += 1

    return palette


@st.cache_data(show_spinner=False)
def get_link_combinations(df_pred, parasite):
    '''
    The predicted interactions of one parasite at the level they can be compared between
    hosts: the orthology group of the parasite protein against the orthology group of
    the host protein.
    '''
    edges = df_pred[df_pred['taxid1_label'] == parasite]
    links = edges.groupby(['group1', 'group2']).agg(
        hosts=('host', lambda h: frozenset(h)),
        interactions=('target', 'size'),
        # tuples, not lists: streamlit hashes a dataframe through pandas, which cannot
        # factorize unhashable values
        parasite_proteins=('source_name', lambda n: tuple(sorted(set(n)))),
        host_proteins=('target_name', lambda n: tuple(sorted(set(str(x).upper() for x in n)))),
        weight=('weight', 'max')).reset_index()
    links['n_hosts'] = links['hosts'].map(len)
    links['combination'] = links['hosts'].map(combination_label)

    return links


def group_label(proteins, group):
    '''
    The name an orthology group is drawn under, on either axis of the matrix: the names
    its proteins carry, since the group id names nothing to read.
    '''
    if not proteins:
        return group

    named = []
    for protein in proteins[:SYMBOLS_IN_LABEL]:
        # the first name goes in whatever its length, so a group is never drawn under an
        # ellipsis alone
        if named and len(', '.join(named + [protein])) > LABEL_CHARS:
            break
        named.append(protein)
    label = ', '.join(named)

    return f'{label}…' if len(named) < len(proteins) else label


def label_margin(labels):
    '''
    How much room the labels along the foot of the matrix need, which plotly does not
    work out for itself: the ticks are forced (see generate_link_matrix), so a label
    longer than the margin is drawn over the edge of the figure rather than dropped.
    '''
    longest = max((len(str(label)) for label in labels), default=0)

    return int(min(220, 40 + 4.3 * longest))


def order_combinations(links):
    '''
    The combinations in the order the figures read them: everything shared first, then
    the smaller sets, and within a size the biggest first, so a figure is read from
    "carried over everywhere" on the left to "one host only" on the right.
    '''
    counts = links.groupby('combination').agg(links=('group2', 'size'),
                                              n_hosts=('n_hosts', 'first'))

    return list(counts.sort_values(['n_hosts', 'links'], ascending=[False, False],
                                   kind='stable').index)


@st.cache_data(show_spinner=False)
def generate_combination_bars(links, all_hosts, config):
    '''
    How many of the transferred interactions carried over to which hosts: a bar per set
    of hosts, over a matrix saying which hosts the set is.
    '''
    order = order_combinations(links)
    counts = links.groupby('combination').size().reindex(order)
    members = links.drop_duplicates('combination').set_index('combination')['hosts']
    palette = combination_palette(links, all_hosts, config)
    colors = [palette[c] for c in order]

    figure = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.05,
                           row_heights=[0.7, 0.3])
    figure.add_trace(go.Bar(x=order, y=counts.values, marker_color=colors,
                            text=counts.values, textposition='outside',
                            hovertemplate='%{x}<br>%{y} orthology-group links<extra></extra>'),
                     row=1, col=1)

    # the matrix below the bars: a filled dot where the set includes the host, a line
    # joining the hosts of a set
    for i, combination in enumerate(order):
        included = [h for h in all_hosts if h in members[combination]]
        if len(included) > 1:
            figure.add_trace(go.Scatter(x=[combination, combination],
                                        y=[included[0], included[-1]], mode='lines',
                                        line=dict(color=colors[i], width=2),
                                        hoverinfo='skip', showlegend=False), row=2, col=1)
        for host in all_hosts:
            inside = host in members[combination]
            figure.add_trace(go.Scatter(x=[combination], y=[host], mode='markers',
                                        marker=dict(size=13, color=colors[i] if inside
                                                    else ABSENT_COLOR),
                                        hoverinfo='skip', showlegend=False), row=2, col=1)

    # explicit margins: plotly cuts off the host names on the rows and the count over the
    # tallest bar
    figure.update_layout(height=460, plot_bgcolor='white', showlegend=False,
                         margin=dict(l=150, r=10, t=45, b=20), bargap=0.4)
    web_utils.count_ticks(figure, counts.max(), axis='y', title_text='orthology-group links',
                showgrid=True, gridcolor='#f0f0f0', row=1, col=1)
    figure.update_yaxes(title_text=None, categoryorder='array', categoryarray=all_hosts[::-1],
                        showgrid=False, row=2, col=1)
    figure.update_xaxes(showticklabels=False, row=1, col=1)
    # the dots below name the set; names run into each other from three hosts on
    figure.update_xaxes(title_text=None, showticklabels=False, showgrid=False, row=2, col=1)

    return figure


@st.cache_data(show_spinner=False)
def generate_link_matrix(links, all_hosts, config):
    '''
    A tile matrix of every transferred interaction: parasite proteins run across, host
    orthology groups run down, and a tile's colour says which hosts received that link.
    '''
    squares = links.copy()
    squares['family'] = [group_label(p, g) for p, g in
                         zip(squares['host_proteins'], squares['group2'])]
    squares['parasite family'] = [group_label(p, g) for p, g in
                                  zip(squares['parasite_proteins'], squares['group1'])]
    # the same gene symbol can occur in separate orthology groups, so the group id is the
    # categorical coordinate on both axes
    squares['family_id'] = squares['group2']
    squares['parasite_id'] = squares['group1']
    squares['host proteins'] = squares['host_proteins'].map(', '.join)
    squares['parasite proteins'] = squares['parasite_proteins'].map(', '.join)
    squares['orthology groups'] = squares['group1'] + ' → ' + squares['group2']
    squares['predicted in'] = squares['hosts'].map(
        lambda hosts: ', '.join(sorted(hosts)))

    host_order = {
        host_label(taxid, config): index
        for index, taxid in enumerate(config['hosts'])
    }
    rows = squares.groupby('family_id').agg(
        family=('family', 'first'),
        hosts=('hosts', lambda values: frozenset().union(*values)),
        links=('family_id', 'size'),
    )
    rows['n_hosts'] = rows['hosts'].map(len)
    rows['host_order'] = rows['hosts'].map(
        lambda hosts: tuple(sorted(host_order[host] for host in hosts)))
    rows = rows.sort_values(['n_hosts', 'host_order', 'links'],
                            ascending=[False, True, False], kind='stable')
    columns = squares.groupby('parasite_id').agg(
        family=('parasite family', 'first'),
        hosts=('n_hosts', 'max'),
        links=('parasite_id', 'size'),
    ).sort_values(['hosts', 'links'], ascending=False, kind='stable')

    order = order_combinations(links)
    palette = combination_palette(links, all_hosts, config)

    figure = px.scatter(squares, x='parasite_id', y='family_id', color='combination',
                        color_discrete_map=palette,
                        # plotly express flips category_orders on a y axis, so the most
                        # shared families come out on top
                        category_orders={'parasite_id': list(columns.index),
                                         'family_id': list(rows.index),
                                         'combination': order},
                        hover_data={'family_id': False, 'parasite_id': False,
                                    'family': True, 'parasite family': False,
                                    'host proteins': True, 'parasite proteins': True,
                                    'orthology groups': True,
                                    'interactions': True, 'combination': False,
                                    'predicted in': True})
    figure.update_traces(marker=dict(symbol='square', size=14,
                                     line=dict(color='white', width=1.5)))
    figure.update_layout(height=max(420, 17 * len(rows) + 260), plot_bgcolor='white',
                         # the names are long enough that plotly cuts them off if the
                         # margins are left to it
                         margin=dict(l=190, r=10, t=10, b=label_margin(columns['family'])),
                         legend=dict(orientation='h', yanchor='bottom', y=1.01, x=0),
                         xaxis_title='parasite protein family',
                         yaxis_title='host protein family',
                         legend_title_text='predicted in host(s)')
    # tickmode='array': plotly thins the labels of an axis this long
    figure.update_xaxes(tickangle=-60, showgrid=True, gridcolor='#eef1f4',
                         gridwidth=1, zeroline=False, tickfont_size=9,
                         tickmode='array',
                         tickvals=list(columns.index), ticktext=list(columns['family']))
    figure.update_yaxes(showgrid=True, gridcolor='#eef1f4', gridwidth=1, zeroline=False,
                         tickfont_size=9, tickmode='array',
                         tickvals=list(rows.index), ticktext=list(rows['family']))

    return figure


@st.cache_data(show_spinner=False)
def count_available_proteins(data_dir, config, parasite_taxid, hosts_taxids):
    '''
    How many host proteins the pipeline had to work with for this parasite in each host,
    which is the number the rest of this page has to be read against.

    The pool is the half of each host the parasite's niche reaches, and not the whole of
    what came through the filters. pipeline/main.py filters the hosts once for every
    parasite, while homology.get_links hands each parasite only its niche's half at
    transfer time, so the whole pool is not what an extracellular parasite was drawn
    from -- for the extracellular parasites of this page it is around three times what
    they were offered, and by a different factor in each host.
    '''
    infected = web_utils.infected_tissue_proteins(data_dir, config, parasite_taxid)
    niche = web_utils.parasite_niche(config, parasite_taxid)

    available = {}
    for host, taxids in hosts_taxids.items():
        proteins = web_utils.filtered_pool(data_dir, taxids, niche=niche)
        available[host] = (len(proteins), len(proteins & infected))

    return available


@st.cache_data(show_spinner=False)
def get_host_coverage(df_pred, parasite, data_dir, config, hosts_taxids):
    '''
    The size of what was predicted in each host, beside the size of what could have
    been.
    '''
    edges = df_pred[df_pred['taxid1_label'] == parasite]
    parasite_taxid = edges['taxid1'].iloc[0]
    available = count_available_proteins(data_dir, config, parasite_taxid, hosts_taxids)

    coverage = edges.groupby('host').agg(**{
        'predicted interactions': ('target', 'size'),
        'parasite proteins': ('source', 'nunique'),
        'host proteins reached': ('target', 'nunique'),
        'host families reached': ('group2', 'nunique')}).reset_index()
    coverage['host proteins this parasite could reach'] = coverage['host'].map(
        lambda h: available[h][0])
    coverage['of those, in a tissue it infects'] = coverage['host'].map(
        lambda h: available[h][1])

    return coverage.rename(columns={'host': 'Host'}).sort_values(
        'predicted interactions', ascending=False, kind='stable')


@st.cache_data(show_spinner=False)
def explain_host_specific(links, all_hosts, data_dir, config, parasite_taxid, hosts_taxids):
    '''
    For every interaction predicted in some hosts but not others, why it is missing from
    the others.
    '''
    orthologs = load_host_orthologs(data_dir)
    if orthologs is None:
        return None

    expressed = web_utils.infected_tissue_proteins(data_dir, config, parasite_taxid)
    # the niche half of the pool, which homology.get_links applies as it transfers and the
    # pool the rest of the page counts against does not
    in_reach = web_utils.filtered_pool(
        data_dir, tuple(taxid for taxids in hosts_taxids.values() for taxid in taxids),
        niche=web_utils.parasite_niche(config, parasite_taxid))

    members, expressed_of, reachable = {}, {}, {}
    for host, taxids in hosts_taxids.items():
        rows = orthologs[orthologs['taxid'].isin(taxids)]
        proteins = rows.groupby('group')['proteins'].apply(
            lambda p: {protein for entry in p for protein in entry.split(',')})
        members[host] = proteins.to_dict()
        expressed_of[host] = {group: proteins_of & expressed
                              for group, proteins_of in members[host].items()}
        reachable[host] = {group: proteins_of & in_reach
                           for group, proteins_of in expressed_of[host].items()}

    missing = []
    for link in links[links['n_hosts'] < len(all_hosts)].itertuples():
        for host in all_hosts:
            if host in link.hosts:
                continue
            family = members[host].get(link.group2, set())
            if not family:
                reason = ABSENT_FAMILY
            elif not expressed_of[host].get(link.group2):
                reason = NOT_EXPRESSED
            elif not reachable[host].get(link.group2):
                reason = OUT_OF_REACH
            else:
                reason = NOT_TRANSFERRED
            missing.append([group_label(link.host_proteins, link.group2),
                            ', '.join(link.parasite_proteins), link.combination, host,
                            len(family), reason])

    return pd.DataFrame(missing, columns=['Host protein family', 'Parasite proteins',
                                          'Predicted in', 'Missing from',
                                          'Proteins of the family in that host', 'Why'])


@st.cache_data(show_spinner=False)
def get_overview(df_pred, config):
    '''
    Every parasite that has more than one host, and how much of it carried over between
    them.
    '''
    overview = []
    for parasite, edges in df_pred.groupby('taxid1_label'):
        hosts = sorted(edges['host'].unique())
        links = get_link_combinations(df_pred, parasite)
        counts = edges.groupby('host').size()
        overview.append([parasite, ', '.join(hosts),
                         ' · '.join(f'{short_name(h)} {counts[h]}' for h in hosts),
                         int((links['n_hosts'] == len(hosts)).sum()),
                         int((links['n_hosts'] == 1).sum())])

    return pd.DataFrame(overview, columns=['Parasite', 'Hosts', 'Predicted interactions',
                                           'Links in every host', 'Links in one host only'])


st.caption('One parasite across the hosts it is predicted against: which interactions '
           'carried over to every host, which are found in a single host, and how the '
           'predictions compare with the host proteins available in each.')
st.markdown("---")

settings, _ = st.columns([1, 1])
with settings:
    score = st.slider('Confidence score', MIN_SCORE, MAX_SCORE, DEFAULT_SCORE,
                      help='Interactions predicted below this confidence are left out, as '
                           'on the other pages. These are the smallest interactomes of the '
                           'set, so this page starts at the bottom of the range: raising it '
                           'leaves fewer parasites with more than one host.')

df_pred = get_multi_host_predictions(data_dir, config, score)

if df_pred.empty:
    st.text('No parasite is predicted against more than one host at this confidence')
else:
    st.subheader('Parasites with more than one host')
    st.caption('Parasites predicted against more than one host, with the number of '
               'interactions predicted in each. The last two columns count interactions as '
               'pairs of orthology groups, the level at which the hosts can be compared: a '
               'link present in every host carried over, a link present in one host only did '
               'not.')
    st.dataframe(get_overview(df_pred, config), width='stretch', hide_index=True)

    parasites = sorted(df_pred['taxid1_label'].unique())
    with st.columns(3)[1]:
        parasite = st.selectbox('Select a parasite to compare its hosts', parasites,
                                index=None, placeholder='<select>')

    if parasite is not None:
        edges = df_pred[df_pred['taxid1_label'] == parasite]
        all_hosts = sorted(edges['host'].unique())
        hosts_taxids = {host: tuple(sorted(set(rows['taxid2'])))
                        for host, rows in edges.groupby('host')}
        links = get_link_combinations(df_pred, parasite)

        # the comparison first, since a host carrying more interactions may only have been
        # annotated more deeply
        st.subheader('Predicted interactions relative to the available host proteins')
        st.caption('Interactions predicted in each host beside the pool of host proteins '
                   'they were drawn from: those passing the tissue and DeepLoc filters and '
                   'lying where this parasite\'s niche can reach them — the surface of the '
                   'cell for a parasite that stays outside it, the cytosol and nucleus as '
                   'well for one that gets in — and those among them annotated to a tissue '
                   'this parasite infects. The second pool is what the predictions were '
                   'drawn from, so a host with more predicted interactions than another may '
                   'simply have more of it.')
        st.dataframe(get_host_coverage(df_pred, parasite, data_dir, config, hosts_taxids),
                     width='stretch', hide_index=True)

        body_column, sets_column = st.columns([1, 1], gap='large')
        with body_column:
            body_figure.show_body_figure(
                config, data_dir, edges[edges['weight'] >= score],
                tuple(taxid for host in all_hosts for taxid in hosts_taxids[host]),
                shared_color_scale=True, title_as_subheader=True)
        with sets_column:
            st.subheader('Shared and host-specific interactions')
            st.caption(f'Interactions of {parasite} predicted in the set of hosts named by '
                       'the matrix below the bars. Interactions are counted as pairs of '
                       'orthology groups rather than pairs of proteins, since the '
                       'orthologous proteins of two hosts are distinct proteins.')
            st.plotly_chart(generate_combination_bars(links, all_hosts, config),
                            width='stretch')

        st.subheader('Interactions per parasite protein family and host protein family')
        st.caption('One tile per transferred interaction: families of parasite proteins on '
                   'the x axis, families of host proteins on the y axis, coloured by the '
                   'host or host set that received the interaction. Both axes are distinct '
                   'orthology groups, labelled with the proteins they hold and with the '
                   'group ID and the full membership on hover, ordered so shared '
                   'interactions gather in the upper left.')
        st.plotly_chart(generate_link_matrix(links, all_hosts, config), width='stretch')

        reasons = explain_host_specific(links, all_hosts, data_dir, config,
                                        edges['taxid1'].iloc[0], hosts_taxids)
        if reasons is not None and not reasons.empty:
            counts = reasons['Why'].value_counts()
            # `none` rather than `0`: a nought in bold reads as a value the page failed to
            # fill in
            absent = counts.get(ABSENT_FAMILY, 0) or 'none'
            said = [f'**Of the {len(reasons)} host–link combinations missing from a host, '
                    f'{absent} are missing because that host has no protein of the family '
                    'at all** — the only reason of the four that is about the host\'s own '
                    'biology rather than about what is known of it.']
            if counts.get(NOT_EXPRESSED):
                said.append(f'{counts[NOT_EXPRESSED]} are a family the host does have, with no '
                            'protein of it annotated to a tissue this parasite infects, so the '
                            'transfer was never attempted. TISSUES cannot tell a protein that '
                            'is not expressed there from one nobody has looked for, and it '
                            'covers human far more deeply than pig, so these are the limit of '
                            'the annotation at least as much as the reach of the parasite.')
            if counts.get(OUT_OF_REACH):
                said.append(f'{counts[OUT_OF_REACH]} are a family the host does express in an '
                            'infected tissue, whose proteins DeepLoc puts where this parasite '
                            'cannot reach them — inside the cell, for a parasite that stays '
                            'outside it.')
            if counts.get(NOT_TRANSFERRED):
                said.append(f'{counts[NOT_TRANSFERRED]} are a family the host has, expresses '
                            'and the parasite reaches, whose link is not recorded under this '
                            'pair of groups. A host protein belonging to both groups of a link '
                            'is the usual reason: the transfer is then written with the two '
                            'groups the other way round, and the host has the interaction '
                            'after all.')
            st.markdown(' '.join(said))
            with st.expander(f'The {len(reasons)} missing links, one row each'):
                st.dataframe(reasons, width='stretch', hide_index=True)
        elif reasons is None:
            st.caption('Run `python scripts/build_host_orthologs.py` to add the check of '
                       'whether a host missing an interaction has a protein of the family '
                       'at all.')

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()
