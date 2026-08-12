import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils
import web_utils
import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from css import style

st.set_page_config(layout="wide", page_title="OrthoHPI 2.0", menu_items={})
st.session_state.data_dir = 'data'
st.session_state.config_file = 'config.yml'
style.load_css()

web_utils.show_pages_menu('Home')

# Read dataset
config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()

# fallback for a parasite without a `group` in the config
UNKNOWN_GROUP = 'Unclassified'
UNKNOWN_COLOR = '#999999'
# A host with two parasites still has to fit two labels under its column, so the columns
# are sized as if every host had at least this many parasites. Human has 35 against the
# two of pig, and strictly proportional columns leave pig a strip its labels overrun.
MIN_COLUMN = 4


def short_name(parasite):
    '''`Plasmodium vivax` as `P. vivax`, the same abbreviation the circos labels its
    arcs with, so the same parasite reads the same way on every page.'''
    return f'{parasite[0]}. {parasite.split(" ")[1]}'


@st.cache_data(show_spinner=False)
def get_overview_predictions(data_dir, config):
    '''
    Every predicted interaction, labelled with the host group it was predicted against
    and with the taxonomic group of its parasite. A parasite infecting more than one host
    appears once per host, since the point of this page is comparing those.

    The hosts come from web_utils.get_host_groups, so the columns of the figures are
    exactly the host options the other two pages offer -- rat and mouse pooled as Rodent.

    :param str data_dir: directory holding predictions.parquet
    :param dict config: parsed configuration
    :return: one row per predicted interaction, with host, parasite, group and weight
    '''
    predictions = web_utils.load_predictions(data_dir)
    groups = {p['label']: p.get('group', UNKNOWN_GROUP) for p in config['parasites'].values()}
    order = {g: i for i, g in enumerate(config.get('parasite_groups', {}))}

    frames = []
    for host, taxids in web_utils.get_host_groups(config, predictions).items():
        frame = predictions.loc[predictions['taxid2'].isin(taxids), ['taxid1_label', 'weight']]
        frames.append(frame.assign(host=host))
    df = pd.concat(frames, ignore_index=True)

    df['group'] = df['taxid1_label'].map(lambda p: groups.get(p, UNKNOWN_GROUP))
    # what puts the parasites of a host in the order of the circos: taxonomic group as
    # the configuration declares it, then name, with an unclassified parasite last
    df['group_rank'] = df['group'].map(lambda g: order.get(g, len(order)))
    df['name'] = df['taxid1_label'].map(short_name)

    return df


def host_columns(df):
    '''
    The skeleton both figures are drawn on: one column per host, as wide as the number of
    parasites infecting it, sharing a y axis so a bar or a box can be compared straight
    across the hosts rather than only within one.

    :param df: overview predictions, as get_overview_predictions builds them
    :return: (figure, [(host, its predictions, its parasites in axis order), ...])
    '''
    hosts = []
    for host in df['host'].unique():
        host_df = df[df['host'] == host]
        parasites = host_df[['taxid1_label', 'name', 'group_rank']].drop_duplicates()
        parasites = parasites.sort_values(by=['group_rank', 'taxid1_label'], kind='stable')
        hosts.append((host, host_df, parasites['name'].tolist()))

    widths = [max(len(names), MIN_COLUMN) for _, _, names in hosts]
    figure = make_subplots(rows=1, cols=len(hosts), shared_yaxes=True,
                           column_widths=[w / sum(widths) for w in widths],
                           subplot_titles=[host for host, _, _ in hosts],
                           horizontal_spacing=0.02)

    return figure, hosts


def style_host_columns(figure, y_title):
    '''The layout the two figures share: the parasite names under each column, the
    quantity named once down the left, and the taxonomic groups as the legend.'''
    figure.update_layout(height=470, plot_bgcolor='white',
                         margin=dict(l=0, r=0, t=95, b=10),
                         legend=dict(orientation='h', yanchor='bottom', y=1.14, x=0,
                                     title_text='', font=dict(size=11)))
    figure.update_xaxes(tickangle=-60, showgrid=False, tickfont=dict(size=11))
    figure.update_yaxes(showgrid=True, gridcolor='#f0f0f0', zerolinecolor='#e0e0e0')
    figure.update_yaxes(title_text=y_title, row=1, col=1)

    return figure


@st.cache_data(show_spinner=False)
def generate_interactions_per_parasite(df, palette):
    '''
    How many interactions are predicted for each parasite, in one column per host.

    The bars are coloured by taxonomic group rather than one colour per parasite: the
    colour then carries something -- the same something the circos, the heatmap strips
    and the dot matrix are coloured by -- instead of only telling neighbouring bars apart,
    which their position already does.
    '''
    figure, hosts = host_columns(df)
    labelled = set()
    for column, (host, host_df, names) in enumerate(hosts, start=1):
        counts = host_df.groupby(['name', 'group'], observed=True).size().reset_index(name='count')
        for group, rows in counts.groupby('group', observed=True):
            figure.add_trace(
                go.Bar(x=rows['name'], y=rows['count'], name=group,
                       marker_color=palette.get(group, UNKNOWN_COLOR),
                       # one legend entry per group however many hosts it turns up in,
                       # and clicking it hides that group in every column at once
                       legendgroup=group, showlegend=group not in labelled,
                       hovertemplate='%{x}<br>%{y} predicted interactions'
                                     f'<extra>{host}</extra>'),
                row=1, col=column)
            labelled.add(group)
        figure.update_xaxes(categoryorder='array', categoryarray=names, row=1, col=column)

    return style_host_columns(figure, 'predicted interactions')


@st.cache_data(show_spinner=False)
def generate_confidence_per_parasite(df, palette):
    '''
    The spread of the confidence score of each parasite's predicted interactions, in the
    same columns and colours as the counts above, so the two figures are read together:
    a parasite with many interactions and a low box has many weakly supported ones.

    The outlying points are left off. Every interaction is a point, and a few hundred of
    them beside each box hide the boxes themselves.
    '''
    figure, hosts = host_columns(df)
    labelled = set()
    for column, (host, host_df, names) in enumerate(hosts, start=1):
        for group, rows in host_df.groupby('group', observed=True):
            figure.add_trace(
                go.Box(x=rows['name'], y=rows['weight'], name=group,
                       marker_color=palette.get(group, UNKNOWN_COLOR),
                       line_width=1, boxpoints=False,
                       legendgroup=group, showlegend=group not in labelled,
                       hovertemplate='%{x}<br>score %{y}' f'<extra>{host}</extra>'),
                row=1, col=column)
            labelled.add(group)
        figure.update_xaxes(categoryorder='array', categoryarray=names, row=1, col=column)

    return style_host_columns(figure, 'confidence score')


st.markdown("<h1 style='text-align: center; color: #023858;'>OrthoHPI 2.0</h1>", unsafe_allow_html=True)
st.markdown("<h3 style='text-align: center; color: #2b8cbe;'>Orthology Prediction of Host-Parasite PPI</h3>", unsafe_allow_html=True)

st.markdown("<h3 style='text-align: center; color: black;'>Every host and every parasite</h3>", unsafe_allow_html=True)
st.caption('Protein-protein interactions between a parasite and its host, predicted by '
           'transferring the interactions of their orthologues. This page is everything that '
           'was predicted, host by host; **Parasites in a host** compares the parasites of one '
           'host with each other, and **Predicted Host-parasite PPIs** opens the network of a '
           'single host and parasite.')
st.markdown("---")

overview = get_overview_predictions(data_dir, config)
parasite_palette = config.get('parasite_groups', {})

st.subheader("Number of Predicted Interactions per Parasite")
st.caption('Every predicted interaction of every parasite, grouped by the host it was '
           'predicted against and coloured by the taxonomic group of the parasite. Nothing is '
           'filtered here -- no confidence threshold and no restriction to the tissues a '
           'parasite infects -- so these counts are the whole prediction and are higher than '
           'the ones the other two pages show, which do restrict to the infected tissues.')
st.plotly_chart(generate_interactions_per_parasite(overview, parasite_palette), width='stretch')

st.subheader("Confidence of the Predicted Interactions per Parasite")
st.caption('How well supported those interactions are: the box spans the middle half of a '
           "parasite's confidence scores, the line across it is the median. The score comes "
           'from the evidence behind the orthologous interaction the prediction was '
           'transferred from, so a high box is a parasite whose predictions rest on '
           'experimental and database evidence rather than on a single weak source.')
st.plotly_chart(generate_confidence_per_parasite(overview, parasite_palette), width='stretch')

# the DeepLoc figures -- where the host proteins that survived the localisation filter sit
# in the cell -- belong here, once data/deeploc is summarised into the app's data directory

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()
