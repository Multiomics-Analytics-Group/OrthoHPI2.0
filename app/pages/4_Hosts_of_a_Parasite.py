import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import utils
import web_utils
import streamlit as st
import textwrap
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import body_figure
from css import style

st.set_page_config(layout="wide", page_title="OrthoHPI 2.0", menu_items={})
style.load_css()
web_utils.show_header('Hosts of a parasite')

# Read dataset
config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()


# The confidence the figures start at, and the range the slider spans, the same as the
# other two pages. Worth knowing while reading it: the parasites with more than one host
# are the small interactomes of the set, and once the predictions are restricted to the
# tissues a parasite infects the usual 0.7 leaves two parasites on the page and one of them
# a single host, so this page opens at the bottom of the range and is raised rather than
# lowered
MIN_SCORE, MAX_SCORE, DEFAULT_SCORE = 0.4, 0.9, 0.4
# One colour per host, fixed here and not read from config['hosts'], which is where the
# rest of the app takes them from. This page is the only one that draws the hosts against
# each other rather than one at a time, and the config colours do not survive that: human
# is a dark grey there, which is the colour of a host a set does not include (ABSENT_COLOR)
# and close to the colour of a link shared by every host. These are the colourblind-safe
# Okabe-Ito set the config already uses for the parasite groups, with the blues left out --
# blue says how far a link carried over on this page, and nothing else may use it.
# A host without an entry falls back to its config colour (see host_palette).
HOST_COLORS = {'Homo sapiens (human)': '#D55E00', 'Sus scrofa (pig)': '#CC79A7',
               'Mus musculus (mouse)': '#009E73', 'Rattus norvegicus (rat)': '#E69F00'}
# colour of a link predicted against every host of the parasite, and the shades of it a
# link predicted against some but not all of them is drawn in -- only reachable with three
# hosts or more, where there is more than one way of being partly shared, so they are
# separated within the blue rather than sharing one colour. A link found in one host only
# is drawn in that host's own colour, from HOST_COLORS
SHARED_COLOR = '#08519c'
PARTIAL_COLORS = ['#6baed6', '#9ecae1', '#4292c6', '#c6dbef']
# how far toward white each host after the first sharing a colour with another is mixed,
# for hosts falling back to their config colour: the config gives rat and mouse the same
# one, since they are one clade and every other page reads them as one host
SPECIES_TINT = 0.45
# what a host is drawn in on the combination matrix when the combination does not include
# it. No host is drawn in a grey, so an empty dot cannot be read as a host colour
ABSENT_COLOR = '#e0e0e0'
# most gene symbols written into the label of an orthology group before the rest are left
# to the hover. A group is a family, not a gene, and the families here run to 18 proteins
SYMBOLS_IN_LABEL = 3
# a GO term is only worth comparing between the two sets of host proteins if it is neither
# nearly empty nor nearly everything. The size is read off every annotated protein of the
# hosts being compared -- two or more of them, so 36,000 to 54,000 proteins, and a term of
# 500 covers fewer than 2 in every 100. Read off the pool the pipeline's filters passed
# instead, the same terms are small enough that a ceiling on the share of it lets the
# broadest of them back in: a quarter of a 2,420-protein pool is 605 proteins, which
# admits `Regulation of primary metabolic process` and fills the figure with terms both
# sets carry alike. This is a count and not an enrichment, and it keeps its own bounds
GO_MIN_PROTEINS, GO_MAX_PROTEINS = 10, 500
# GO terms drawn in the comparison of what the shared and the host-specific proteins do,
# and the width their names are broken over lines at beside the axis -- the same width the
# network page wraps them at
GO_TOP_N = 12
GO_LABEL_WRAP_WIDTH = 42
# how many host proteins of a set have to carry a GO term for it to be drawn at all
GO_MIN_IN_SET = 2
# what the two sets of host proteins are called wherever they are counted or drawn
SHARED_SET, SPECIFIC_SET = 'shared by every host', 'host-specific'
def short_name(parasite):
    '''`Trichinella spiralis` as `T. spiralis`, the abbreviation the other pages use.'''
    parts = parasite.split(' ')

    return f'{parasite[0]}. {parts[1]}' if len(parts) > 1 else parasite


@st.cache_data(show_spinner=False)
def load_host_orthologs(data_dir):
    '''
    Which proteins of each host belong to each of the orthology groups the predictions
    reach, written by scripts/build_host_orthologs.py.

    It answers the question the predictions cannot: a group with no predicted interaction
    in pig is missing from predictions.parquet whether pig has no protein in it at all or
    whether pig's proteins were filtered out before the transfer, and those two are not
    the same finding. A data directory built before that script simply leaves the section
    that reads this out.

    :param str data_dir: directory holding host_orthologs.parquet
    :return: dataframe of group, taxid, n_proteins, proteins; or None if not built
    '''
    input_file = os.path.join(data_dir, 'host_orthologs.parquet')
    if not os.path.exists(input_file):
        return None

    return utils.read_parquet_file(input_file=input_file)


def host_label(host, config):
    '''The configured species name of a host taxid.'''
    return config['hosts'][int(host)]['label']


def host_color(host_label, config):
    '''The colour this page fixes for the host (HOST_COLORS), falling back to the one the
    config gives it -- a host added to the config and not to HOST_COLORS still gets drawn
    in the colour the rest of the app knows it by.'''
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

    Every prediction is kept whatever tissue the host protein is annotated to. The tissue
    filter of the pipeline keeps a host protein expressed in a tissue *any* parasite
    infects, so a host protein can carry a predicted interaction with a parasite that
    never reaches the tissue it was kept for; the figure at the foot of the page is where
    that is read, and restricting the whole page to it would empty two of the comparisons.

    :param float score: interactions predicted below this confidence are left out
    :return: predictions of the multi-host parasites, with a `host` column
    '''
    predictions = web_utils.load_predictions(data_dir)
    predictions = predictions[predictions['weight'] >= score].copy()
    predictions['host'] = predictions['taxid2'].map(lambda t: host_label(t, config))

    hosts = predictions.groupby('taxid1')['host'].nunique()

    return predictions[predictions['taxid1'].isin(hosts[hosts > 1].index)]


def combination_label(hosts):
    '''The name of a set of hosts: the hosts abbreviated and joined, in a fixed order so
    that the same set is always the same label and the same column of the figures.'''
    return ' + '.join(short_name(host) for host in sorted(hosts))


def tint(color, amount):
    '''
    Mixes a colour toward white, the same way the network page lightens a species colour
    into a node fill.

    :param str color: '#rrggbb'
    :param float amount: 0 leaves the colour alone, 1 turns it white
    :return: the mixed colour, or the colour unchanged if it is not a hex triplet
    '''
    color = str(color)
    if not color.startswith('#') or len(color) != 7:
        return color
    channels = [int(color[i:i + 2], 16) for i in (1, 3, 5)]

    return '#%02x%02x%02x' % tuple(round(c + (255 - c) * amount) for c in channels)


def host_palette(all_hosts, config):
    '''
    One colour per host of the parasite. A repeated configured colour is tinted for later
    hosts so each species remains distinguishable.

    :return: {host label: colour}
    '''
    sharing = {}
    for host in all_hosts:
        sharing.setdefault(host_color(host, config), []).append(host)

    return {host: tint(color, min(SPECIES_TINT * i, 0.7))
            for color, hosts in sharing.items() for i, host in enumerate(hosts)}


def combination_palette(links, all_hosts, config):
    '''
    The colour each set of hosts is drawn in, over both figures, so that a set is the same
    colour wherever it is read.

    The colouring says how far a link carried over and not which host it is: everything
    shared is one colour, a link in one host only is that host's own, and the sets in
    between -- which only exist for a parasite with three hosts -- take shades of the
    shared colour, one each, since three of them in one blue cannot be told apart in a
    legend.

    :return: {combination label: colour}
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
    hosts: the orthology group of the parasite protein against the orthology group of the
    host protein.

    A human protein and a pig protein are different proteins and cannot be intersected,
    which is why counting host proteins per host says nothing about whether the same
    interaction was predicted twice. The pair of orthology groups can be: it is what the
    pipeline transferred the link at (homology.get_links reads the COG links between two
    groups), so two hosts carrying the same pair carry the same transferred interaction,
    and a pair in one host only is an interaction that did not carry over.

    Two things are worth keeping in mind reading it. A group is a protein family and not
    a gene -- KOG0994 is LAMA1, LAMB1 and LAMB2 in human -- so a shared pair means the
    same family was reached, not necessarily the same gene. And the host proteins of the
    pair are named per host, so the label of a row is the union of what the hosts call it.

    :param df_pred: multi-host predictions
    :param str parasite: label of the parasite to read
    :return: one row per (parasite group, host group) pair, with the hosts it was
             predicted in, the proteins on both sides, and the parasite protein names
    '''
    edges = df_pred[df_pred['taxid1_label'] == parasite]
    links = edges.groupby(['group1', 'group2']).agg(
        hosts=('host', lambda h: frozenset(h)),
        interactions=('target', 'size'),
        # tuples and not lists: the frame is passed into the cached figure functions, and
        # streamlit hashes a dataframe through pandas, which cannot factorize a column of
        # unhashable values and falls back to pickling the whole frame
        parasite_proteins=('source_name', lambda n: tuple(sorted(set(n)))),
        host_proteins=('target_name', lambda n: tuple(sorted(set(str(x).upper() for x in n)))),
        weight=('weight', 'max')).reset_index()
    links['n_hosts'] = links['hosts'].map(len)
    links['combination'] = links['hosts'].map(combination_label)

    return links


def group_label(host_proteins, group):
    '''The name a host orthology group is drawn under: the gene symbols its host proteins
    carry, since the group id names nothing to read. The symbols are the union over the
    hosts, cut where the axis stops earning the space -- the whole of it is in the hover
    -- and the id is kept as the label of a group whose proteins have no symbol.'''
    if not host_proteins:
        return group
    label = ', '.join(host_proteins[:SYMBOLS_IN_LABEL])

    return f'{label}…' if len(host_proteins) > SYMBOLS_IN_LABEL else label


def order_combinations(links):
    '''The combinations in the order the figures read them: everything shared first, then
    the smaller sets, and within a size the biggest first, so a figure is read from
    "carried over everywhere" on the left to "one host only" on the right.'''
    counts = links.groupby('combination').agg(links=('group2', 'size'),
                                              n_hosts=('n_hosts', 'first'))

    return list(counts.sort_values(['n_hosts', 'links'], ascending=[False, False],
                                   kind='stable').index)


@st.cache_data(show_spinner=False)
def generate_combination_bars(links, all_hosts, config):
    '''
    How many of the transferred interactions carried over to which hosts: a bar per set of
    hosts, over a matrix saying which hosts the set is.

    A Venn diagram is the usual figure for this and does not survive a third host; the
    bars do, and they also put the sets in an order -- shared first, one host only last --
    that a Venn has no way of showing.

    :param links: the link combinations of one parasite
    :param list all_hosts: every host of the parasite, in the order of the matrix rows
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

    # the matrix below the bars: a filled dot where the set includes the host, an empty
    # one where it does not, and a line joining the hosts of a set that spans several
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

    # the left margin is explicit rather than left to plotly, which cuts off the host names
    # labelling the rows of the matrix, and the top one leaves room for the count over the
    # tallest bar
    figure.update_layout(height=460, plot_bgcolor='white', showlegend=False,
                         margin=dict(l=150, r=10, t=45, b=20), bargap=0.4)
    web_utils.count_ticks(figure, counts.max(), axis='y', title_text='orthology-group links',
                showgrid=True, gridcolor='#f0f0f0', row=1, col=1)
    figure.update_yaxes(title_text=None, categoryorder='array', categoryarray=all_hosts[::-1],
                        showgrid=False, row=2, col=1)
    figure.update_xaxes(showticklabels=False, row=1, col=1)
    # the dots below name the set, and a name spelled out under them as well runs into its
    # neighbour as soon as a parasite has three hosts
    figure.update_xaxes(title_text=None, showticklabels=False, showgrid=False, row=2, col=1)

    return figure


@st.cache_data(show_spinner=False)
def generate_link_matrix(links, all_hosts, config):
    '''
    Every transferred interaction of the parasite as one square: the parasite protein
    across, the host protein family down, the colour saying which hosts it was predicted
    in.

    Two networks side by side is the other way of drawing this and is the wrong one: the
    reader has to hold one of them in their head to find what the other is missing, and
    the host proteins have different names in each. Here the difference is the figure --
    a row that changes colour across is a family one host reaches and another does not.

    The rows are ordered by how many hosts reach them, then by host in configured order;
    the block of shared interactions gathers in the top left corner and host-specific
    rows stay together instead of mixing Rat and Mouse families.
    '''
    squares = links.explode('parasite_proteins').rename(
        columns={'parasite_proteins': 'parasite protein'})
    squares['family'] = [group_label(p, g) for p, g in
                         zip(squares['host_proteins'], squares['group2'])]
    squares['host proteins'] = squares['host_proteins'].map(', '.join)
    squares['orthology groups'] = squares['group1'] + ' → ' + squares['group2']

    host_order = {
        host_label(taxid, config): index
        for index, taxid in enumerate(config['hosts'])
    }
    rows = squares.groupby('family').agg(
        hosts=('hosts', lambda values: frozenset().union(*values)),
        links=('family', 'size'),
    )
    rows['n_hosts'] = rows['hosts'].map(len)
    rows['host_order'] = rows['hosts'].map(
        lambda hosts: tuple(sorted(host_order[host] for host in hosts)))
    rows = rows.sort_values(['n_hosts', 'host_order', 'links'],
                            ascending=[False, True, False], kind='stable')
    columns = (squares.groupby('parasite protein')
                      .agg(hosts=('n_hosts', 'max'), links=('parasite protein', 'size'))
                      .sort_values(['hosts', 'links'], ascending=False, kind='stable'))

    order = order_combinations(links)
    palette = combination_palette(links, all_hosts, config)

    figure = px.scatter(squares, x='parasite protein', y='family', color='combination',
                        color_discrete_map=palette,
                        # plotly express flips category_orders on a y axis, so the most
                        # shared families first here puts them in the top rows
                        category_orders={'parasite protein': list(columns.index),
                                         'family': list(rows.index),
                                         'combination': order},
                        hover_data={'family': False, 'host proteins': True,
                                    'orthology groups': True, 'interactions': True,
                                    'combination': True})
    figure.update_traces(marker=dict(symbol='square', size=9, line=dict(width=0)))
    figure.update_layout(height=max(420, 17 * len(rows) + 260), plot_bgcolor='white',
                         # the gene symbols down the side and the parasite proteins along the
                         # foot are long enough that plotly cuts them off if left to itself
                         margin=dict(l=190, r=10, t=10, b=120), legend_title_text='',
                         legend=dict(orientation='h', yanchor='bottom', y=1.01, x=0),
                         xaxis_title='protein of the parasite',
                         yaxis_title='host protein family')
    figure.update_xaxes(tickangle=-60, showgrid=True, gridcolor='#f7f7f7')
    figure.update_yaxes(showgrid=True, gridcolor='#f7f7f7')

    return figure


@st.cache_data(show_spinner=False)
def count_available_proteins(data_dir, config, parasite_taxid, hosts_taxids):
    '''
    How many host proteins the pipeline had to work with in each host, which is the
    number the rest of this page has to be read against.

    Two counts per host: the proteins that came through the secretome, tissue and DeepLoc
    filters at all, and the ones among them annotated to a tissue this parasite is known
    to infect. Both come from tissues_cell_types.parquet, which pipeline/main.py writes
    for exactly the proteins that survived the filters.

    :param dict hosts_taxids: {host label: tuple of taxids as str}
    :return: {host label: (proteins after the filters, of those in an infected tissue)}
    '''
    tissues = utils.read_parquet_file(input_file=f'{data_dir}/tissues_cell_types.parquet')
    infected = web_utils.infected_tissue_proteins(data_dir, config, parasite_taxid)

    available = {}
    for host, taxids in hosts_taxids.items():
        proteins = set(tissues.loc[tissues['Gene'].str.split('.').str[0].isin(taxids), 'Gene'])
        available[host] = (len(proteins), len(proteins & infected))

    return available


@st.cache_data(show_spinner=False)
def get_host_coverage(df_pred, parasite, data_dir, config, hosts_taxids):
    '''
    The size of what was predicted in each host, beside the size of what could have been.

    The point of the table is that the first set of numbers is not a result on its own.
    Human carries two to three times the interactions of pig for the same parasite, and
    the last two columns are why: the human proteome is annotated far more deeply, so far
    more human proteins reach the transfer step at all.

    The two 'available' columns are the pool the pipeline drew from, and the second of them
    is a superset of the proteins reached: both are restricted to the tissues this parasite
    infects, so the two can be read against each other.
    '''
    edges = df_pred[df_pred['taxid1_label'] == parasite]
    parasite_taxid = edges['taxid1'].iloc[0]
    available = count_available_proteins(data_dir, config, parasite_taxid, hosts_taxids)

    coverage = edges.groupby('host').agg(**{
        'predicted interactions': ('target', 'size'),
        'parasite proteins': ('source', 'nunique'),
        'host proteins reached': ('target', 'nunique'),
        'host families reached': ('group2', 'nunique')}).reset_index()
    coverage['host proteins available after the filters'] = coverage['host'].map(
        lambda h: available[h][0])
    coverage['available in a tissue this parasite infects'] = coverage['host'].map(
        lambda h: available[h][1])

    return coverage.rename(columns={'host': 'Host'}).sort_values(
        'predicted interactions', ascending=False, kind='stable')


@st.cache_data(show_spinner=False)
def explain_host_specific(links, all_hosts, data_dir, config, parasite_taxid, hosts_taxids):
    '''
    For every interaction predicted in some hosts but not others, why it is missing from
    the others. There are three ways it can be, and they mean different things:

    * the host has no protein in that orthology group at all -- the family is absent from
      the host, and the parasite has nothing to interact with. This is the only one of the
      three that is about the host's biology.
    * the host has proteins in the group, but none of them came through the pipeline's
      filters into a tissue the parasite infects, so the transfer was never attempted.
    * the host has such a protein and the link was still not transferred. One run of the
      pipeline cannot produce this -- the proteins that came through the filters are
      exactly the ones the orthology groups were matched from, so one of them in the group
      would have carried the link -- but a data directory whose parquet files were written
      by different runs can.

    On the data as it stands every host-specific link of every parasite is the second,
    which is worth saying plainly: the differences between the hosts on this page are
    differences in how deeply the hosts are annotated, not in what the parasite can reach.
    The check is kept rather than hard-coded because that stops being true as soon as a
    host or a parasite is added.

    :return: (dataframe of one row per missing link, None if host_orthologs is not built)
    '''
    orthologs = load_host_orthologs(data_dir)
    if orthologs is None:
        return None

    expressed = web_utils.infected_tissue_proteins(data_dir, config, parasite_taxid)

    members, reachable = {}, {}
    for host, taxids in hosts_taxids.items():
        rows = orthologs[orthologs['taxid'].isin(taxids)]
        proteins = rows.groupby('group')['proteins'].apply(
            lambda p: {protein for entry in p for protein in entry.split(',')})
        members[host] = proteins.to_dict()
        reachable[host] = {group: proteins_of & expressed
                           for group, proteins_of in members[host].items()}

    missing = []
    for link in links[links['n_hosts'] < len(all_hosts)].itertuples():
        for host in all_hosts:
            if host in link.hosts:
                continue
            family = members[host].get(link.group2, set())
            if not family:
                reason = 'no protein of the family in this host'
            elif not reachable[host].get(link.group2):
                reason = 'family present, but no protein of it reaches the infected tissues'
            else:
                reason = 'protein available, link not transferred'
            missing.append([group_label(link.host_proteins, link.group2),
                            ', '.join(link.parasite_proteins), link.combination, host,
                            len(family), reason])

    return pd.DataFrame(missing, columns=['Host protein family', 'Parasite proteins',
                                          'Predicted in', 'Missing from',
                                          'Proteins of the family in that host', 'Why'])


@st.cache_data(show_spinner=False)
def get_go_comparison(data_dir, shared_proteins, specific_proteins, taxids):
    '''
    What the host proteins of the shared interactions do, beside what the host-specific
    ones do.

    A count and not an enrichment: a host-specific set here runs to three proteins, and a
    Fisher test on three proteins reports whatever the smallest set happens to contain.
    The terms are filtered to GO_MIN_PROTEINS..GO_MAX_PROTEINS of the annotated proteins
    of the hosts being compared, which is what keeps `Cellular process` out of a figure
    meant to separate two sets. Those are this figure's own bounds and not the ones the
    network page's enrichment uses, which are a share of a much smaller background.

    :param tuple shared_proteins: host proteins carrying an interaction found in every host
    :param tuple specific_proteins: host proteins carrying one found in some hosts only
    :param tuple taxids: host taxids the GO annotations are read for
    :return: dataframe of term, set, proteins; or None if nothing passes the filters
    '''
    species = [int(t) for t in taxids]
    go_df = utils.read_parquet_file(input_file=f'{data_dir}/gos.parquet',
                                    filters=[('taxid', 'in', species)])
    go_df = go_df[go_df['taxid'].isin(species)]

    sizes = go_df.groupby('description')['#string_protein_id'].nunique()
    general = set(sizes[(sizes < GO_MIN_PROTEINS) | (sizes > GO_MAX_PROTEINS)].index)

    counted = []
    for name, proteins in ((SHARED_SET, shared_proteins),
                           (SPECIFIC_SET, specific_proteins)):
        terms = go_df[go_df['#string_protein_id'].isin(proteins)]
        terms = terms[~terms['description'].isin(general)]
        counts = terms.groupby('description')['#string_protein_id'].nunique()
        counts = counts[counts >= GO_MIN_IN_SET]
        # the size of the set the count came out of, which is what makes a count of two
        # comparable between a set of sixty proteins and a set of three
        counted.append(pd.DataFrame({'term': counts.index, 'set': name,
                                     'proteins': counts.values,
                                     'set_size': len(proteins)}))

    comparison = pd.concat(counted, ignore_index=True)

    return comparison if not comparison.empty else None


@st.cache_data(show_spinner=False)
def generate_go_bars(comparison):
    '''
    The GO terms of the two sets of host proteins as paired bars.

    Only the terms both sets carry are drawn, and they are ranked by how different the
    share of each set carrying them is. Ranking by the count instead let the larger set
    choose the terms: the shared set is two to five times the host-specific one for most
    parasites, so the twelve biggest counts were twelve terms of the shared set, drawn
    against nothing and read as a comparison. Ranking by the difference alone is no better,
    since the largest difference is a term one set does not carry at all.

    :param comparison: term, set, proteins and set_size, as get_go_comparison builds it
    :return: the figure, or None if no term is carried by both sets
    '''
    counts = comparison.pivot_table(index='term', columns='set', values='proteins',
                                    fill_value=0)
    for name in (SHARED_SET, SPECIFIC_SET):
        if name not in counts.columns:
            return None

    both = counts[(counts[SHARED_SET] > 0) & (counts[SPECIFIC_SET] > 0)]
    if both.empty:
        return None

    sizes = comparison.drop_duplicates('set').set_index('set')['set_size']
    difference = (both[SHARED_SET] / sizes[SHARED_SET]
                  - both[SPECIFIC_SET] / sizes[SPECIFIC_SET]).abs()
    terms = list(difference.sort_values(ascending=False, kind='stable').head(GO_TOP_N).index)
    view = comparison[comparison['term'].isin(terms)].copy()
    # the names of GO terms are sentences, and one of them beside an axis is as wide as the
    # figure, so they are broken over lines as the network page breaks them
    wrapped = {term: '<br>'.join(textwrap.wrap(str(term), width=GO_LABEL_WRAP_WIDTH))
               for term in terms}
    view['term'] = view['term'].map(wrapped)

    figure = px.bar(view, x='proteins', y='term', color='set', orientation='h',
                    barmode='group',
                    color_discrete_map={SHARED_SET: SHARED_COLOR,
                                        SPECIFIC_SET: PARTIAL_COLORS[0]},
                    category_orders={'term': [wrapped[t] for t in terms]})
    figure.update_layout(height=max(360, 46 * len(terms) + 140), plot_bgcolor='white',
                         margin=dict(l=330, r=10, t=10, b=40), legend_title_text='',
                         legend=dict(orientation='h', yanchor='bottom', y=1.01, x=0),
                         xaxis_title='host proteins annotated to the term',
                         yaxis_title=None, bargap=0.3)
    web_utils.count_ticks(figure, view['proteins'].max(), showgrid=True, gridcolor='#f0f0f0')

    return figure


@st.cache_data(show_spinner=False)
def get_overview(df_pred, config):
    '''
    Every parasite that has more than one host, and how much of it carried over between
    them. Drawn before anything is selected, since which parasites are even comparable is
    the first thing to know -- there are seven of them as species, four with rat and mouse
    read as one host.
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

# The confidence threshold decides which multi-host predictions there are to compare.
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

        body_column, coverage_column = st.columns([1, 1], gap='large')
        with body_column:
            body_figure.show_body_figure(
                config, data_dir, edges[edges['weight'] >= score],
                tuple(taxid for host in all_hosts for taxid in hosts_taxids[host]),
                shared_color_scale=True, title_as_subheader=True)
        with coverage_column:
            st.subheader('Predicted interactions relative to the available host proteins')
            st.caption('Interactions predicted in each host beside the pool of host proteins '
                       'available: those passing the secretome, tissue and DeepLoc filters, '
                       'and those among them annotated to a tissue this parasite infects. The '
                       'second pool is what the predictions were drawn from, so a host with '
                       'more predicted interactions than another may simply have more of it.')
            st.dataframe(
                get_host_coverage(df_pred, parasite, data_dir, config, hosts_taxids),
                width='stretch', hide_index=True)

        # side by side: the bars are three columns wide whatever the parasite, which is a
        # figure that does not need the width of the page, and the matrix beside them says
        # what the sets of the bars are made of
        sets, matrix = st.columns([2, 3], gap='large')
        with sets:
            st.subheader('Shared and host-specific interactions')
            st.caption(f'Interactions of {parasite} predicted in the set of hosts named by '
                       'the matrix below the bars. Interactions are counted as pairs of '
                       'orthology groups rather than pairs of proteins, since the '
                       'orthologous proteins of two hosts are distinct proteins. See the '
                       'caveat below the figures before interpreting a host-specific bar as '
                       'host specificity.')
            st.plotly_chart(generate_combination_bars(links, all_hosts, config),
                            width='stretch')

        with matrix:
            st.subheader('Interactions per parasite protein and host protein family')
            st.caption('One square per predicted interaction: parasite proteins on the x '
                       'axis, families of host proteins on the y axis, coloured by the hosts '
                       'the interaction was predicted in. Rows are orthology groups rather '
                       'than genes, named on hover, and are ordered so that the shared '
                       'interactions gather in the upper left.')
            st.plotly_chart(generate_link_matrix(links, all_hosts, config),
                            width='stretch')

        reasons = explain_host_specific(links, all_hosts, data_dir, config,
                                        edges['taxid1'].iloc[0], hosts_taxids)
        if reasons is not None and not reasons.empty:
            counts = reasons['Why'].value_counts()
            absent = counts.get('no protein of the family in this host', 0)
            st.markdown(f'**Of the {len(reasons)} host–link combinations missing from a host, '
                        f'{absent} are missing because that host has no protein of the family '
                        'at all.** Every other one is a family the host does have, whose '
                        'proteins were filtered out before the transfer was attempted — so '
                        'it is a gap in the host\'s annotation and not a protein the parasite '
                        'cannot reach.')
            with st.expander(f'The {len(reasons)} missing links, one row each'):
                st.dataframe(reasons, width='stretch', hide_index=True)
        elif reasons is None:
            st.caption('Run `python scripts/build_host_orthologs.py` to add the check of '
                       'whether a host missing an interaction has a protein of the family '
                       'at all.')

        shared = tuple(edges.loc[edges['group2'].isin(
            links.loc[links['n_hosts'] == len(all_hosts), 'group2']), 'target'].unique())
        specific = tuple(edges.loc[edges['group2'].isin(
            links.loc[links['n_hosts'] == 1, 'group2']), 'target'].unique())
        taxids = tuple(sorted(set(edges['taxid2'])))
        comparison = get_go_comparison(data_dir, shared, specific, taxids)
        go_bars = generate_go_bars(comparison) if comparison is not None else None
        with st.container():
            st.subheader('Gene Ontology terms of shared and host-specific proteins')
            if go_bars is not None:
                st.caption('Gene Ontology terms carried by the host proteins of the '
                           'interactions found in every host and by those of the '
                           'interactions found in a single host, for the terms both sets '
                           'carry and ranked by how different the share of each set is. '
                           'These are counts of proteins and not an enrichment. Terms '
                           f'carried by fewer than {GO_MIN_PROTEINS} or more than '
                           f'{GO_MAX_PROTEINS} of the hosts\' annotated proteins are '
                           'omitted, as neither separates the two sets.')
                st.plotly_chart(go_bars, width='stretch')
            elif comparison is not None:
                st.caption('No term is carried by both the shared and the host-specific '
                           'proteins of this parasite, so there is nothing to compare them '
                           'on. The host-specific set is the smaller of the two and runs to '
                           'a handful of proteins.')
            else:
                st.caption('None of the host proteins of this parasite carries a Gene '
                           'Ontology term in the size range the figure selects on.')

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()
