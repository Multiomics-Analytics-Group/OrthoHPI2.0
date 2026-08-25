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
from css import style

st.set_page_config(layout="wide", page_title="OrthoHPI 2.0", menu_items={})
style.load_css()
web_utils.show_pages_menu('Hosts of a parasite')

# Read dataset
config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()


# The confidence the figures start at, and the range the slider spans, the same as the
# other two pages. Worth knowing while reading it: the parasites with more than one host
# are the small interactomes of the set, and the default leaves some of them three
# interactions in a host -- lower the slider to see the comparison at full size
MIN_SCORE, MAX_SCORE, DEFAULT_SCORE = 0.4, 0.9, 0.7
# One colour per host, fixed here and not read from config['hosts'], which is where the
# rest of the app takes them from. This page is the only one that draws the hosts against
# each other rather than one at a time, and the config colours do not survive that: human
# is a dark grey there, which is the colour of a host a set does not include (ABSENT_COLOR)
# and close to the colour of a link shared by every host, and rat and mouse are the one
# green of their clade, which is two of the columns here. These are the colourblind-safe
# Okabe-Ito set the config already uses for the parasite groups, with the blues left out --
# blue says how far a link carried over on this page, and nothing else may use it.
# A host without an entry falls back to its config colour (see host_palette).
HOST_COLORS = {'Homo sapiens': '#D55E00', 'Sus scrofa': '#CC79A7',
               'Mus musculus': '#009E73', 'Rattus norvegicus': '#E69F00',
               'Rodent': '#009E73'}
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
# nearly empty nor nearly everything; the same bounds utils.calculate_enrichment selects on
GO_MIN_PROTEINS, GO_MAX_PROTEINS = 10, 500
# GO terms drawn in the comparison of what the shared and the host-specific proteins do,
# and the width their names are broken over lines at beside the axis -- the same width the
# network page wraps them at
GO_TOP_N = 12
GO_LABEL_WRAP_WIDTH = 42
# how many host proteins of a set have to carry a GO term for it to be drawn at all
GO_MIN_IN_SET = 2


def short_name(parasite):
    '''`Trichinella spiralis` as `T. spiralis`, the abbreviation the other pages label
    their parasites with. Host names go through it too, so `Rodent` -- a host group that
    is a clade and not a binomial -- is returned unchanged rather than split.'''
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


def host_label(host, config, pooled):
    '''The name a host taxid is read under on this page: its own species name, or the
    group the config collapses it into (rat and mouse -> Rodent) when the pooling box is
    ticked, which is how every other page reads it.'''
    host = config['hosts'][int(host)]

    return host.get('group', host['label']) if pooled else host['label']


def host_color(host_label, config, pooled):
    '''The colour this page fixes for the host (HOST_COLORS), falling back to the one the
    config gives it -- a host added to the config and not to HOST_COLORS still gets drawn,
    in the colour the rest of the app knows it by. Grouped hosts carry the same colour in
    the config, so a group has one whichever of its members is found.'''
    if host_label in HOST_COLORS:
        return HOST_COLORS[host_label]

    for taxid, host in config['hosts'].items():
        if (host.get('group', host['label']) if pooled else host['label']) == host_label:
            return host['color']

    return SHARED_COLOR


@st.cache_data(show_spinner=False)
def get_multi_host_predictions(data_dir, config, pooled, score=MIN_SCORE):
    '''
    Every prediction of the parasites that are predicted against more than one host,
    labelled with the host it was predicted against.

    Whether rat and mouse count as one host or two is the whole population of this page:
    five of the seven parasites with more than one host are rodent parasites, and pooling
    them the way the rest of the app does leaves four parasites to compare. Both readings
    are offered, and the parasites that qualify are recounted under the one chosen.

    Every prediction is kept whatever tissue the host protein is annotated to. The tissue
    filter of the pipeline keeps a host protein expressed in a tissue *any* parasite
    infects, so a host protein can carry a predicted interaction with a parasite that
    never reaches the tissue it was kept for; the figure at the foot of the page is where
    that is read, and restricting the whole page to it would empty two of the comparisons.

    :param bool pooled: read the hosts as the config groups them rather than as species
    :param float score: interactions predicted below this confidence are left out
    :return: predictions of the multi-host parasites, with a `host` column
    '''
    predictions = web_utils.load_predictions(data_dir)
    predictions = predictions[predictions['weight'] >= score].copy()
    predictions['host'] = predictions['taxid2'].map(lambda t: host_label(t, config, pooled))

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


def host_palette(all_hosts, config, pooled):
    '''
    One colour per host of the parasite, which the config cannot always give: rat and
    mouse carry the same colour there because they are one clade and every other page
    reads them as one host. Read as species -- which is how this page reads them by
    default, since five of the seven parasites with more than one host are rodent
    parasites -- that is two hosts in one colour and two columns that cannot be told apart.

    The first host of a colour keeps it and the others are mixed toward white, so the
    clade is still one colour and the species within it are still separate.

    :return: {host label: colour}
    '''
    sharing = {}
    for host in all_hosts:
        sharing.setdefault(host_color(host, config, pooled), []).append(host)

    return {host: tint(color, min(SPECIES_TINT * i, 0.7))
            for color, hosts in sharing.items() for i, host in enumerate(hosts)}


def combination_palette(links, all_hosts, config, pooled):
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
    hosts_palette = host_palette(all_hosts, config, pooled)
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
def generate_combination_bars(links, all_hosts, config, pooled):
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
    palette = combination_palette(links, all_hosts, config, pooled)
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
    figure.update_yaxes(title_text='orthology-group links', showgrid=True,
                        gridcolor='#f0f0f0', row=1, col=1)
    figure.update_yaxes(title_text=None, categoryorder='array', categoryarray=all_hosts[::-1],
                        showgrid=False, row=2, col=1)
    figure.update_xaxes(showticklabels=False, row=1, col=1)
    # the dots below name the set, and a name spelled out under them as well runs into its
    # neighbour as soon as a parasite has three hosts
    figure.update_xaxes(title_text=None, showticklabels=False, showgrid=False, row=2, col=1)

    return figure


@st.cache_data(show_spinner=False)
def generate_link_matrix(links, all_hosts, config, pooled):
    '''
    Every transferred interaction of the parasite as one square: the parasite protein
    across, the host protein family down, the colour saying which hosts it was predicted
    in.

    Two networks side by side is the other way of drawing this and is the wrong one: the
    reader has to hold one of them in their head to find what the other is missing, and
    the host proteins have different names in each. Here the difference is the figure --
    a row that changes colour across is a family one host reaches and another does not.

    The rows are ordered by how many hosts reach them and the columns by how many rows
    they carry, so the block of shared interactions gathers in the top left corner and
    what only one host has falls away from it.
    '''
    squares = links.explode('parasite_proteins').rename(
        columns={'parasite_proteins': 'parasite protein'})
    squares['family'] = [group_label(p, g) for p, g in
                         zip(squares['host_proteins'], squares['group2'])]
    squares['host proteins'] = squares['host_proteins'].map(', '.join)
    squares['orthology groups'] = squares['group1'] + ' → ' + squares['group2']

    rows = (squares.groupby('family')
                   .agg(hosts=('n_hosts', 'max'), links=('family', 'size'))
                   .sort_values(['hosts', 'links'], ascending=False, kind='stable'))
    columns = (squares.groupby('parasite protein')
                      .agg(hosts=('n_hosts', 'max'), links=('parasite protein', 'size'))
                      .sort_values(['hosts', 'links'], ascending=False, kind='stable'))

    order = order_combinations(links)
    palette = combination_palette(links, all_hosts, config, pooled)

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
    figure.update_traces(marker=dict(symbol='square', size=11, line=dict(width=0)))
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
def get_infected_tissue_proteins(data_dir, config, parasite_taxid):
    '''
    The host proteins annotated to a tissue this parasite is known to infect, over every
    host at once.

    The pipeline's tissue filter is not parasite-specific -- it keeps a host protein
    expressed in a tissue *any* parasite of the config infects -- so a host protein can
    carry a predicted interaction with a parasite that never reaches the tissue it was
    kept for. This is the set that says which ones do not.

    :return: set of STRING protein ids
    '''
    tissues = utils.read_parquet_file(input_file=f'{data_dir}/tissues_cell_types.parquet')
    mapped = config['tissues']
    infected = {mapped[t].lower() for t in config['parasites'][int(parasite_taxid)]['tissues']}

    return set(tissues.loc[tissues['Tissue'].isin(infected), 'Gene'])


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
    infected = get_infected_tissue_proteins(data_dir, config, parasite_taxid)

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

    The two 'available' columns are the pool the pipeline drew from and are not a superset
    of the proteins reached: the tissue filter is not parasite-specific, so a host protein
    kept for a tissue another parasite infects can still carry an interaction here. That is
    why what was reached is counted twice -- in any tissue, and in a tissue this parasite
    infects -- rather than compared against the pool.
    '''
    edges = df_pred[df_pred['taxid1_label'] == parasite]
    parasite_taxid = edges['taxid1'].iloc[0]
    available = count_available_proteins(data_dir, config, parasite_taxid, hosts_taxids)
    infected = get_infected_tissue_proteins(data_dir, config, parasite_taxid)

    coverage = edges.groupby('host').agg(**{
        'predicted interactions': ('target', 'size'),
        'parasite proteins': ('source', 'nunique'),
        'host proteins reached': ('target', 'nunique'),
        'host families reached': ('group2', 'nunique')}).reset_index()
    reached_here = edges[edges['target'].isin(infected)].groupby('host')['target'].nunique()
    coverage['reached in a tissue this parasite infects'] = coverage['host'].map(
        reached_here).fillna(0).astype(int)
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

    expressed = get_infected_tissue_proteins(data_dir, config, parasite_taxid)

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
    The terms are filtered to the size range utils.calculate_enrichment selects on, which
    is what keeps `Cellular process` out of a figure meant to separate two sets.

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
    for name, proteins in (('shared by every host', shared_proteins),
                           ('host-specific', specific_proteins)):
        terms = go_df[go_df['#string_protein_id'].isin(proteins)]
        terms = terms[~terms['description'].isin(general)]
        counts = terms.groupby('description')['#string_protein_id'].nunique()
        counts = counts[counts >= GO_MIN_IN_SET]
        counted.append(pd.DataFrame({'term': counts.index, 'set': name,
                                     'proteins': counts.values}))

    comparison = pd.concat(counted, ignore_index=True)

    return comparison if not comparison.empty else None


@st.cache_data(show_spinner=False)
def generate_go_bars(comparison):
    '''The GO terms of the two sets of host proteins as paired bars, the terms ranked by
    how many proteins carry them in either set. A term drawn against one set only is one
    the other set has fewer than GO_MIN_IN_SET proteins for, not one it has none for.'''
    ranked = comparison.groupby('term')['proteins'].max().sort_values(ascending=False,
                                                                     kind='stable')
    terms = list(ranked.head(GO_TOP_N).index)
    view = comparison[comparison['term'].isin(terms)].copy()
    # the names of GO terms are sentences, and one of them beside an axis is as wide as the
    # figure, so they are broken over lines as the network page breaks them
    wrapped = {term: '<br>'.join(textwrap.wrap(str(term), width=GO_LABEL_WRAP_WIDTH))
               for term in terms}
    view['term'] = view['term'].map(wrapped)

    figure = px.bar(view, x='proteins', y='term', color='set', orientation='h',
                    barmode='group',
                    color_discrete_map={'shared by every host': SHARED_COLOR,
                                        'host-specific': PARTIAL_COLORS[0]},
                    category_orders={'term': [wrapped[t] for t in terms]})
    figure.update_layout(height=max(360, 46 * len(terms) + 140), plot_bgcolor='white',
                         margin=dict(l=330, r=10, t=10, b=40), legend_title_text='',
                         legend=dict(orientation='h', yanchor='bottom', y=1.01, x=0),
                         xaxis_title='host proteins annotated to the term',
                         yaxis_title=None, bargap=0.3)
    figure.update_xaxes(showgrid=True, gridcolor='#f0f0f0')

    return figure


@st.cache_data(show_spinner=False)
def count_interactions_per_host_tissue(df_pred, parasite, data_dir, config):
    '''
    Where in each host the predicted interactions can take place: the interactions whose
    host protein is annotated to a tissue the parasite is known to infect, counted per
    host and tissue.

    This is the one place the tissue restriction is applied, and it is why the rest of the
    page does not apply it: a host whose proteins are annotated to none of the tissues the
    parasite infects -- rat, for two of the rodent parasites -- disappears from it, and a
    comparison with one host left in it is not a comparison.
    '''
    edges = df_pred[df_pred['taxid1_label'] == parasite]
    tissues = utils.read_parquet_file(input_file=f'{data_dir}/tissues_cell_types.parquet')
    tissues = tissues.rename({'Gene': 'target'}, axis=1)[['target', 'Tissue']].drop_duplicates()

    mapped = config['tissues']
    infected = {mapped[t].lower() for t in config['parasites'][int(edges['taxid1'].iloc[0])]['tissues']}

    aux = edges[['host', 'source', 'target']].drop_duplicates()
    aux = pd.merge(aux, tissues, on='target')
    aux = aux[aux['Tissue'].isin(infected)]
    if aux.empty:
        return None

    return aux.groupby(['host', 'Tissue']).size().rename('interactions').reset_index()


@st.cache_data(show_spinner=False)
def generate_host_tissue_dots(per_tissue, all_hosts, config, pooled):
    '''A dot wherever a host carries a predicted interaction with a protein expressed in a
    tissue the parasite infects, sized by how many. A gap in a row is a host whose proteins
    of that tissue are not annotated, which on this data is what most of the gaps are.'''
    dots = per_tissue.copy()
    tissues = list(dots.groupby('Tissue')['interactions'].sum().sort_values(
        ascending=False, kind='stable').index)
    palette = host_palette(all_hosts, config, pooled)

    figure = px.scatter(dots, x='host', y='Tissue', color='host', size='interactions',
                        size_max=26, color_discrete_map=palette,
                        category_orders={'host': all_hosts, 'Tissue': tissues},
                        hover_data={'host': True, 'Tissue': True, 'interactions': True})
    figure.update_traces(marker=dict(sizemin=5, line=dict(width=0)))
    # a column per host and nothing else, so the figure carries its own width rather than
    # being stretched to the page: two hosts spread over a metre of it are two columns with
    # a field between them. The tissue names down the side need the left margin spelled out
    figure.update_layout(height=max(320, 34 * len(tissues) + 200), plot_bgcolor='white',
                         width=170 * len(all_hosts) + 200, showlegend=False,
                         margin=dict(l=170, r=20, t=10, b=40),
                         xaxis_title=None, yaxis_title='tissue the parasite infects')
    figure.update_xaxes(showgrid=True, gridcolor='#f0f0f0')
    figure.update_yaxes(showgrid=True, gridcolor='#f0f0f0')

    return figure


@st.cache_data(show_spinner=False)
def get_overview(df_pred, config, pooled):
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


st.markdown("<h1 style='text-align: center; color: #023858;'>OrthoHPI 2.0</h1>", unsafe_allow_html=True)
st.markdown("<h3 style='text-align: center; color: #2b8cbe;'>Orthology Prediction of Host-Parasite PPIs</h3>", unsafe_allow_html=True)

st.markdown("<h3 style='text-align: center; color: black;'>One parasite, several hosts</h3>", unsafe_allow_html=True)
st.markdown("---")

# the two settings that decide what there is to compare sit above everything, since both
# of them change which parasites are even on the page
settings, _ = st.columns([1, 1])
with settings:
    pooled = st.checkbox('Read rat and mouse as one host (Rodent), as the other pages do',
                         value=False,
                         help='Five of the seven parasites with more than one host are '
                              'rodent parasites. Pooling rat and mouse the way the rest of '
                              'the app does leaves four parasites to compare.')
    score = st.slider('Confidence score', MIN_SCORE, MAX_SCORE, DEFAULT_SCORE,
                      help='Interactions predicted below this confidence are left out, as '
                           'on the other pages. These are the smallest interactomes of the '
                           'set, so lowering it is often what makes the comparison big '
                           'enough to read.')

df_pred = get_multi_host_predictions(data_dir, config, pooled, score)

if df_pred.empty:
    st.text('No parasite is predicted against more than one host at this confidence')
else:
    st.subheader('The Parasites With More Than One Host')
    st.caption('Every parasite the pipeline predicts against more than one host, with how '
               'many interactions were predicted in each. The last two columns count them '
               'at the level the hosts can be compared at -- the orthology group of the '
               'parasite protein against the orthology group of the host protein, which is '
               'what the pipeline transferred the interaction at -- so a link in every host '
               'is one that carried over, and a link in one host only is one that did not.')
    st.dataframe(get_overview(df_pred, config, pooled), width='stretch', hide_index=True)

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

        st.subheader('Which Interactions Carried Over to Which Hosts')
        st.caption(f'Each bar is a set of hosts and counts the interactions of {parasite} '
                   'predicted in exactly those hosts, with the set named by the matrix '
                   'below it. An interaction is counted as a pair of orthology groups '
                   'rather than as a pair of proteins, since a human protein and a pig '
                   'protein are different proteins and cannot be intersected; the pair of '
                   'groups is what the interaction was transferred at, so the same pair in '
                   'two hosts is the same transferred interaction. Read the caveat below '
                   'the figures before reading a host-specific bar as host specificity.')
        st.plotly_chart(generate_combination_bars(links, all_hosts, config, pooled),
                        width='stretch')

        st.subheader('Where the Hosts Differ, Protein by Protein')
        st.caption('One square per predicted interaction: a protein of the parasite across, '
                   'the family of host proteins it reaches down, and the colour saying which '
                   'hosts it was predicted in. A row that is one colour all the way across is '
                   'a family every host reaches; a row that changes colour is where the hosts '
                   'part. The rows are families and not genes -- an orthology group can hold '
                   'several genes of the same host, which the hover names -- and they are '
                   'ordered so that what is shared gathers in the top left.')
        st.plotly_chart(generate_link_matrix(links, all_hosts, config, pooled),
                        width='stretch')

        st.subheader('How Much of the Difference Is the Host, and How Much Is the Annotation')
        st.caption('The size of what was predicted in each host, beside the size of what '
                   'could have been. The last two columns are the pool the pipeline drew '
                   'from: the host proteins that came through its secretome, tissue and '
                   'DeepLoc filters, and the ones among them annotated to a tissue this '
                   'parasite infects. A host with more predicted interactions than another '
                   'has, on this data, more proteins in that pool rather than more of the '
                   'parasite\'s attention. The pool is not a superset of what was reached: '
                   'the tissue filter keeps a host protein expressed in a tissue any '
                   'parasite infects, so what was reached is counted twice instead — in any '
                   'tissue, and in a tissue this parasite infects.')
        st.dataframe(get_host_coverage(df_pred, parasite, data_dir, config, hosts_taxids),
                     width='stretch', hide_index=True)

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
        if comparison is not None:
            st.subheader('What the Shared and the Host-Specific Proteins Do')
            st.caption('The Gene Ontology terms carried by the host proteins of the '
                       'interactions found in every host, beside those carried by the host '
                       'proteins of the interactions found in one host only. These are '
                       'counts of proteins and not an enrichment: the host-specific set '
                       'runs to a handful of proteins, and a Fisher test on a handful '
                       'reports whatever it happens to contain. Terms annotated to fewer '
                       f'than {GO_MIN_PROTEINS} or more than {GO_MAX_PROTEINS} proteins of '
                       'the host are left out, which is what keeps the top of the ontology '
                       'out of the figure.')
            st.plotly_chart(generate_go_bars(comparison), width='stretch')

        per_tissue = count_interactions_per_host_tissue(df_pred, parasite, data_dir, config)
        if per_tissue is not None:
            st.subheader('Where the Interactions Can Take Place in Each Host')
            st.caption('The predicted interactions whose host protein is annotated to a '
                       'tissue this parasite is known to infect, per host and tissue. This '
                       'is the only figure on the page that applies that restriction: the '
                       'tissue filter of the pipeline keeps a host protein expressed in a '
                       'tissue any parasite infects, so elsewhere on the page an interaction '
                       'is counted whatever tissue its host protein was kept for. A host '
                       'missing from a row has no annotated protein there, which is again '
                       'more often the annotation than the biology.')
            st.plotly_chart(generate_host_tissue_dots(per_tissue, all_hosts, config, pooled),
                            width='content')

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()
