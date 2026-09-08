"""
Shades the TISSUES body figure of a host with the number of predicted interactions
reaching each of its organs.

The figures (images/tissues/tissues_<species>.svg) come from tissues.jensenlab.org and draw 21
organs, each element carrying the organ as its `title` attribute. That is a coarser
vocabulary than the 33 lifecycle tissues of the configuration, which is what decides
whether a prediction is shown at all; the organs are read from the same TISSUES download
by scripts/build_figure_tissues.py, so nothing has to be mapped between the two.
"""
import os
import xml.etree.ElementTree as ET
from copy import deepcopy

import pandas as pd
import streamlit as st
try:
    from st_click_detector import click_detector
except ImportError:  # the figure is still drawn, it just cannot be clicked
    click_detector = None

import utils
# the same BTO code -> organ map the annotation was built with, rather than a second copy
# of it here: the two have to name the 21 organs identically or the shading silently
# misses some
from scripts.build_figure_tissues import FIGURE_ORGANS

SVG_NS = 'http://www.w3.org/2000/svg'
ET.register_namespace('', SVG_NS)
ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')

FIGURE_DIR = os.path.join('images', 'tissues')

# The pig says `thyroid` where the other figures say `thyroid gland`; all figures call
# their urinary-bladder region `urine`. Normalize these source labels for annotations and
# user-facing tooltips.
ORGAN_ALIASES = {
    'thyroid': 'thyroid gland',
    'urine': 'urinary bladder',
}

# A blue interaction-count ramp. The first colour is for an organ no interaction reaches,
# so it stays the white the figures are drawn in.
NO_INTERACTIONS_COLOR = '#ffffff'
PALETTE = ['#deebf7', '#9ecae1', '#6baed6', '#3182bd', '#08519c']

# The organs of a clickable figure are anchors, and the component hands back the id of
# the one that was clicked. The organ already standing for the tissue filter carries this
# id instead of its own name, so clicking it a second time clears the filter rather than
# setting what is already set -- and so two clicks on one organ never send the same value
# twice, which is what a click has to do to be noticed.
CLEAR_ORGAN = '__clear__'
# session_state key the click of each host's figure is stored under, per taxid
ORGAN_CLICK_KEY = 'organ_click'
# the outline drawn around the organ the tissue filter is currently set to, and the one a
# clickable organ takes on hover. The organs are filled by interaction count, so the
# selection cannot be another fill without saying something about the count
SELECTED_OUTLINE = '#2b8cbe'
HOVER_OUTLINE = '#7fc0dd'

# the legend is stripped and redrawn below the figure: it is labelled with the confidence
# scores of the TISSUES website, which are not what is being shown here, and the pig
# draws those labels as outlines rather than text, so they cannot simply be rewritten
LEGEND_ID = 'Legend'

# drawing units left between the bottom of the body and the edge of the cropped figure
CROP_MARGIN = 10
# Height of the anatomy row when several hosts are compared. The maps share this height
# where their width permits, keeping the host labels directly above the human figure.
COMPARISON_FIGURE_HEIGHT = 200
# The human source SVG contains frontal and side views on one canvas. The compact
# comparison uses only the frontal view, which covers the organs shown for comparisons.
HUMAN_FRONT_VIEW = (220, 220)
HUMAN_SIDE_VIEW = (600, 125)
HUMAN_VIEW_GAP = 10
COMPACT_HUMAN_FIGURE_HEIGHT = 360

# The tissues a parasite infects are the 33 fine-grained BTO terms of config['tissues'];
# the figure draws 21 coarser organs (20 for rat, which has no gall bladder). Most terms
# are drawn directly by FIGURE_ORGANS. The following lifecycle tissues have a defensible
# coarser organ on the figure: seven parasites are recorded only in `small intestine`,
# which is represented as `intestine`, and would otherwise shade nothing at all.
# (`gastrointestinal tract` is in config['tissues'] but no parasite uses it; it is kept
# here so the map covers the whole vocabulary.)
#
# Five terms still have nothing standing for them:
#   `macrophage`        Leishmania braziliensis, L. infantum, L. major
#   `mouth`             Leishmania braziliensis
#   `nose`              Leishmania braziliensis
#   `placenta`          Toxoplasma gondii
#   `vagina`            Trichomonas vaginalis
# Those terms shade no organ rather than being shaded onto a neighbouring one, which would
# be inventing a location for them. Their interactions still pass the tissue filter and
# appear in the network and the table; they just contribute to no organ on the figure.
# Trichomonas vaginalis is the only parasite whose tissues are all in this list, so it is
# the only one whose figure stays blank -- the others also infect an organ that is drawn.
ORGAN_PARENTS = {
    'BTO:0000142': 'nervous system',      # brain
    'BTO:0001279': 'nervous system',      # spinal cord
    'BTO:0000651': 'intestine',           # small intestine
    'BTO:0000269': 'intestine',           # colon
    'BTO:0001158': 'intestine',           # rectum
    'BTO:0000511': 'intestine',           # gastrointestinal tract
    'BTO:0000122': 'liver',               # bile duct
    'BTO:0000752': 'lymph nodes',         # lymph vessel
    'BTO:0000779': 'intestine',           # mesenteric artery
    'BTO:0001426': 'urinary bladder',     # urethra
}


@st.cache_data(show_spinner=False)
def load_figure_tissues(data_dir, modified_at):
    '''
    Host proteins annotated with the organs of the body figure, written by
    scripts/build_figure_tissues.py. The file is not part of the older snapshot data
    directories, so a missing one only leaves the figure out.

    :param str data_dir: directory holding figure_tissues.parquet
    :param float modified_at: figure_tissues.parquet modification time, used to invalidate
                              cached annotations after a rebuild
    :return: dataframe of Gene and Organ, or None when it has not been built
    '''
    input_file = os.path.join(data_dir, 'figure_tissues.parquet')
    if not os.path.exists(input_file):
        return None

    return utils.read_parquet_file(input_file=input_file)


def get_species(config, taxid):
    '''
    The species name jensenlab uses for a host, taken from the tissue url the
    configuration already holds so that the hosts are only listed in one place. It names
    both the download and the figure, images/tissues/tissues_<species>.svg.

    :param dict config: parsed configuration
    :param taxid: host taxid, as a string or an int
    :return: species name, or None when the host has no tissue annotation
    '''
    host = config['hosts'].get(int(taxid), {})
    url = host.get('tissues_url')
    if url is None:
        return None

    return os.path.basename(url).split('_')[0]


def legend_top(legend):
    '''
    Where the legend of a figure starts, in the coordinates of the drawing, so that the
    space it occupied can be cropped away once it is removed.

    :param legend: the legend group element
    :return: the y it starts at, or None when it cannot be worked out
    '''
    tops = [float(rect.get('y')) for rect in legend.iter(f'{{{SVG_NS}}}rect')
            if rect.get('y') is not None]
    if not tops:
        return None

    # the legends are placed with a plain translate, the only transform the figures use
    offset = 0.0
    transform = legend.get('transform', '')
    if transform.startswith('translate('):
        parts = transform[len('translate('):].rstrip(')').replace(',', ' ').split()
        if len(parts) > 1:
            offset = float(parts[1])

    return min(tops) + offset


def crop_to(root, above):
    '''
    Shortens the viewBox of a figure so it ends just above the given y, which is where
    its legend used to be. Leaving the viewBox alone would keep drawing the empty band
    the legend occupied, pushing the body up its column.

    :param root: the svg root element
    :param above: y the figure should end above, or None to leave the viewBox alone
    '''
    view_box = root.get('viewBox')
    if above is None or view_box is None:
        return

    x, y, width, height = (float(value) for value in view_box.replace(',', ' ').split())
    cropped = above - y - CROP_MARGIN
    if 0 < cropped < height:
        root.set('viewBox', f'{x} {y} {width} {cropped}')


@st.cache_data(show_spinner=False)
def load_figure(species):
    '''
    The body figure of a species, with its legend removed and its fixed size replaced by
    one that follows the column it is drawn in.

    :param str species: species name as jensenlab spells it (human, mouse, rat, pig)
    :return: (svg root element as a string, organs the figure draws), or (None, set())
    '''
    figure_file = os.path.join(FIGURE_DIR, f'tissues_{species}.svg')
    if not os.path.exists(figure_file):
        return None, set()

    root = ET.parse(figure_file).getroot()
    # located before anything is removed, so the walk is not cut short by the tree
    # changing under it
    legends = [(parent, child) for parent in root.iter() for child in parent
               if child.get('id') == LEGEND_ID]
    for parent, legend in legends:
        # The non-human legends overlap lower anatomy in their SVG coordinate space, so
        # shortening their viewBox cuts off the feet, tail, or lower body.
        if species == 'human':
            crop_to(root, above=legend_top(legend))
        parent.remove(legend)

    # the figures are drawn at a fixed pixel size, which would overflow a narrow column
    root.attrib.pop('width', None)
    root.attrib.pop('height', None)
    root.set('style', 'width: 100%; height: auto;')

    organs = {ORGAN_ALIASES.get(element.get('title'), element.get('title'))
              for element in root.iter() if element.get('title')}

    return ET.tostring(root, encoding='unicode'), organs


def infected_organs(config, taxid):
    '''
    The organs of the body figure that a parasite is recorded as infecting.

    The tissue filter that decides which predictions are shown at all keeps a host protein
    expressed in one of those tissues, but TISSUES then annotates that protein to every
    organ it is detected in, most of which the parasite never reaches: a Loa loa protein
    selected for being expressed in skin also comes annotated to the nervous system, and
    the figure drew the brain as the darkest organ on the page. Restricting the shading to
    the organs the parasite actually infects is what keeps the figure about the parasite
    rather than about how broadly its targets happen to be expressed.

    :param dict config: parsed configuration
    :param taxid: parasite taxid, as a string or an int
    :return: set of organ names as the figures label them, possibly empty
    '''
    tissues = config['parasites'].get(int(taxid), {}).get('tissues', [])

    organs = set()
    for code in tissues:
        organ = FIGURE_ORGANS.get(code, ORGAN_PARENTS.get(code))
        if organ is not None:
            organs.add(ORGAN_ALIASES.get(organ, organ))

    return organs


def tissue_organs(config, tissues):
    '''
    The body-figure organs corresponding to explicitly selected lifecycle tissues.

    :param dict config: parsed configuration
    :param iterable tissues: selected, lower-case display names
    :return: set of normalized organ names represented by those tissues
    '''
    selected = {tissue.lower() for tissue in tissues}
    organs = set()
    for code, tissue in config['tissues'].items():
        if tissue.lower() not in selected:
            continue
        organ = FIGURE_ORGANS.get(code, ORGAN_PARENTS.get(code))
        if organ is not None:
            organs.add(ORGAN_ALIASES.get(organ, organ))

    return organs


def organ_tissues(config, organ):
    """
    The lifecycle tissues an organ of the figure stands for, which is what clicking it
    puts in the tissue filter. The inverse of tissue_organs: `liver` selects liver and
    bile duct, `intestine` the whole gut vocabulary.

    :param dict config: parsed configuration
    :param str organ: normalized organ name
    :return: sorted lower-case display names, as the tissue filter offers them
    """
    tissues = []
    for code, tissue in config['tissues'].items():
        drawn = FIGURE_ORGANS.get(code, ORGAN_PARENTS.get(code))
        if drawn is not None and ORGAN_ALIASES.get(drawn, drawn) == organ:
            tissues.append(tissue.lower())

    return sorted(tissues)


def clicked_organ(taxids):
    """
    The organ clicked on one of the host figures since the last run, or None when the
    click is one that has already been acted on.

    Each figure keeps its own component value, which the component leaves in place across
    reruns; an organ is therefore only acted on when the value it is read from changes.
    Every figure's value is marked as seen whichever one carried the click, so a stale
    click on a host that has since been hidden cannot fire later.

    :param iterable taxids: taxids of the hosts whose figures were drawn
    :return: the organ name, CLEAR_ORGAN, or None
    """
    organ = None
    for taxid in taxids:
        key = f'{ORGAN_CLICK_KEY}_{taxid}'
        clicked = st.session_state.get(key)
        if clicked and clicked != st.session_state.get(f'{key}_seen'):
            st.session_state[f'{key}_seen'] = clicked
            organ = organ or clicked

    return organ


def apply_organ_click(config, taxids, options, key):
    """
    Writes the tissues of a clicked organ into the tissue filter, and drops any tissue the
    selected parasite does not reach -- switching parasite otherwise leaves the filter
    holding an option it is no longer offered, which streamlit rejects.

    Called before the filter widget is created, since that is the only point at which its
    value can still be set.

    :param dict config: parsed configuration
    :param iterable taxids: taxids of the hosts whose figures were drawn
    :param iterable options: the tissues the filter offers for this parasite
    :param str key: session_state key of the tissue filter
    """
    organ = clicked_organ(taxids)
    if organ == CLEAR_ORGAN:
        st.session_state[key] = []
    elif organ is not None:
        st.session_state[key] = [tissue for tissue in organ_tissues(config, organ)
                                 if tissue in options]
    elif key in st.session_state:
        st.session_state[key] = [tissue for tissue in st.session_state[key]
                                 if tissue in options]


def link_organs(svg, clickable, selected, config):
    """
    Turns the organs of a figure into links the click component can report, and says in
    the tooltip what clicking one does.

    Every organ the parasite infects is a link, whether or not any interaction reaches it,
    so the filter can be moved straight from one organ to another rather than having to be
    cleared in between.

    :param str svg: the shaded figure
    :param set clickable: organs to link
    :param set selected: organs the tissue filter is currently set to
    :param dict config: parsed configuration
    :return: the figure with its organs wrapped in anchors
    """
    root = ET.fromstring(svg)
    style = ET.Element(f'{{{SVG_NS}}}style')
    style.text = (f'a {{ cursor: pointer; }} '
                  f'a:hover [fill] {{ stroke: {HOVER_OUTLINE}; stroke-width: 2; }}')
    root.insert(0, style)

    # the parents are taken first: the organs are moved inside new elements, and the tree
    # must not be rearranged while it is being walked
    organs = [(parent, child) for parent in root.iter() for child in parent
              if child.get('title')]
    for parent, element in organs:
        organ = ORGAN_ALIASES.get(element.get('title'), element.get('title'))
        if organ not in clickable:
            continue

        if organ in selected:
            for shape in element.iter():
                if shape.get('fill') is not None:
                    shape.set('stroke', SELECTED_OUTLINE)
                    shape.set('stroke-width', '2')
            hint = 'click to clear the tissue filter'
        else:
            hint = f"click to filter on {', '.join(organ_tissues(config, organ))}"

        title = element.find(f'{{{SVG_NS}}}title')
        if title is not None:
            title.text = f'{title.text} -- {hint}'

        link = ET.Element(f'{{{SVG_NS}}}a', {'id': CLEAR_ORGAN if organ in selected else organ})
        parent.insert(list(parent).index(element), link)
        parent.remove(element)
        link.append(element)

    return inline(ET.tostring(root, encoding='unicode'))


def count_interactions(df, figure_tissues):
    '''
    Counts the predicted interactions reaching each organ. An interaction is counted once
    per organ its host protein is annotated to, so the counts do not add up to the size
    of the network -- most host proteins are annotated to a single organ, but a broadly
    expressed one carries as many as twenty.

    The predictions are repeated once per tissue and single-cell cluster of their host
    protein, so they are reduced to one row per interaction first.

    :param df: predictions dataframe, already filtered to what the network shows
    :param figure_tissues: dataframe of Gene and Organ
    :return: {organ: number of interactions}
    '''
    interactions = df.drop_duplicates(subset=['source', 'target'])
    counts = pd.merge(interactions[['target']], figure_tissues,
                      left_on='target', right_on='Gene')

    return counts['Organ'].value_counts().to_dict()


def color_scale(counts):
    '''
    Splits the interaction counts into the shades of the palette, over the range the
    figure actually spans rather than a fixed one, so that a small network is not drawn
    uniformly pale.

    The bins grow geometrically rather than being equal in width. The counts are heavily
    skewed -- `nervous system` collects several times what any other organ does, being
    both the most studied tissue and the one the predictions target most -- and equal-width
    bins put nearly every other organ in the palest one, which reads as if only the brain
    were targeted at all.

    :param dict counts: {organ: number of interactions}
    :return: (list of (inclusive upper bound, colour) from palest, highest count);
             ([], 0) when nothing was counted
    '''
    highest = max(counts.values(), default=0)
    if highest == 0:
        return [], 0

    edges = []
    for i in range(len(PALETTE)):
        # each bin covers the same factor, and never repeats the previous bound: a
        # network whose busiest organ has only a few interactions gets fewer bins
        upper = max(round(highest ** ((i + 1) / len(PALETTE))),
                    edges[-1] + 1 if edges else 1)
        if upper >= highest:
            edges.append(highest)
            break
        edges.append(upper)

    # the colours are spread over the whole palette rather than taken from one end, so
    # the busiest organ is the darkest blue whatever the size of the network
    bounds = []
    for i, upper in enumerate(edges):
        shade = (round(i * (len(PALETTE) - 1) / (len(edges) - 1))
                 if len(edges) > 1 else len(PALETTE) - 1)
        bounds.append((upper, PALETTE[shade]))

    return bounds, highest


def organ_color(count, bounds):
    for upper, color in bounds:
        if count <= upper:
            return color

    return NO_INTERACTIONS_COLOR


def shade_figure(svg, counts, bounds):
    '''
    Colours each organ of the figure by the number of interactions reaching it, and gives
    it a tooltip saying so. The fill sits on the element that carries the organ name, and
    on its children when the organ is drawn as a group of several shapes.

    :param str svg: the figure, as returned by load_figure
    :param dict counts: {organ: number of interactions}
    :param list bounds: colour bins, as returned by color_scale
    :return: the figure with its organs coloured
    '''
    root = ET.fromstring(svg)
    # collected before anything is changed: the tooltips are added as children, and the
    # tree must not grow while it is being walked
    organ_elements = [element for element in root.iter() if element.get('title')]

    for element in organ_elements:
        organ = ORGAN_ALIASES.get(element.get('title'), element.get('title'))
        count = counts.get(organ, 0)
        color = organ_color(count, bounds) if count else NO_INTERACTIONS_COLOR

        # an organ drawn as a group of shapes carries its fill on each of them
        for shape in element.iter():
            if shape is element or shape.get('fill') is not None:
                shape.set('fill', color)

        interactions = 'interaction' if count == 1 else 'interactions'
        title = ET.Element(f'{{{SVG_NS}}}title')
        title.text = f'{organ}: {count} predicted {interactions}'
        element.insert(0, title)

    return inline(ET.tostring(root, encoding='unicode'))


def inline(svg):
    '''
    Puts the whole drawing on one line, which is what it takes for st.markdown to render
    it. The figures were drawn in Illustrator, which wraps long path definitions over
    several indented lines and leaves a blank line between some of them; markdown ends a
    block of raw HTML at the blank line and treats the indented remainder as a code
    block, so the tail of the drawing was printed as text underneath it.

    :param str svg: the drawing
    :return: the same drawing without line breaks
    '''
    return ' '.join(svg.split('\n'))


def bottom_aligned(svg):
    '''Fits an SVG inside the comparison row and aligns its lower edge with the others.'''
    root = ET.fromstring(svg)
    if root.get('id') == 'human':
        return frontal_human_view(root)

    root.set('style', 'width: auto; height: 100%; max-width: 100%;')

    return (f'<div style="height: {COMPARISON_FIGURE_HEIGHT}px; display: flex; '
            f'align-items: flex-end; justify-content: center;">'
            f'{inline(ET.tostring(root, encoding="unicode"))}</div>')


def frontal_human_view(root):
    '''Uses the frontal human anatomy view for the compact multi-host comparison.'''
    _, y, _, height = (float(value) for value in root.get('viewBox').replace(',', ' ').split())
    x, width = HUMAN_FRONT_VIEW
    root.set('viewBox', f'{x} {y} {width} {height}')
    root.set('style', 'width: auto; height: 100%; max-width: 100%;')

    return (f'<div style="height: {COMPARISON_FIGURE_HEIGHT}px; display: flex; '
            f'align-items: flex-end; justify-content: center;">'
            f'{inline(ET.tostring(root, encoding="unicode"))}</div>')


def compact_human_views(svg):
    '''Places cropped frontal and side human views together for the detailed network page.'''
    root = ET.fromstring(svg)
    _, y, _, height = (float(value) for value in root.get('viewBox').replace(',', ' ').split())
    views = []
    for x, width in (HUMAN_FRONT_VIEW, HUMAN_SIDE_VIEW):
        view = deepcopy(root)
        view.set('viewBox', f'{x} {y} {width} {height}')
        view.set('style', 'width: auto; height: 100%; max-width: 100%;')
        views.append(inline(ET.tostring(view, encoding='unicode')))

    return (f'<div style="height: {COMPACT_HUMAN_FIGURE_HEIGHT}px; display: flex; '
            f'align-items: flex-end; justify-content: center; gap: {HUMAN_VIEW_GAP}px;">'
            + ''.join(views) + '</div>')


def legend_html(bounds, compact=False):
    '''The colour bins as a row of swatches, labelled with the counts they stand for.'''
    swatches = []
    previous = 0
    swatch_height = 10 if compact else 14
    label_size = '0.65rem' if compact else '0.7rem'
    for upper, color in bounds:
        label = str(upper) if upper == previous + 1 else f'{previous + 1}-{upper}'
        swatches.append(
            f'<div style="text-align: center; flex: 1;">'
            f'<div style="background: {color}; border: 1px solid #939598; '
            f'height: {swatch_height}px;"></div>'
            f'<div style="font-size: {label_size}; color: #555;">{label}</div></div>')
        previous = upper

    max_width = 'max-width: 240px;' if compact else ''
    return (f'<div style="display: flex; gap: 2px; margin-top: 0.5rem; {max_width}">'
            + ''.join(swatches) + '</div>')


def show_body_figure(config, data_dir, df, taxids, selected_tissues=None,
                     shared_color_scale=False, compact_human=False, title_as_subheader=False,
                     clickable=False):
    '''
    Draws the body figure of each selected host, its organs shaded by the number of
    predicted interactions reaching them. Each selected host species gets its own figure,
    annotated against its own TISSUES data.

    :param dict config: parsed configuration
    :param str data_dir: directory holding figure_tissues.parquet
    :param df: predictions dataframe, already filtered to what the network shows
    :param taxids: taxids of the selected host
    :param iterable selected_tissues: explicitly selected tissue display names, if any
    :param bool shared_color_scale: use one interaction-count scale and legend across all
                                    host figures
    :param bool compact_human: bring the frontal and side human views closer together
    :param bool title_as_subheader: match the surrounding page's section-heading size
    :param bool clickable: draw the organs as links that set the tissue filter. The click
                           is read back by apply_organ_click on the run that follows it
    '''
    figure_tissues_file = os.path.join(data_dir, 'figure_tissues.parquet')
    modified_at = (os.path.getmtime(figure_tissues_file)
                   if os.path.exists(figure_tissues_file) else None)
    figure_tissues = load_figure_tissues(data_dir, modified_at)
    if figure_tissues is None or df.empty:
        return

    drawn = [(taxid, get_species(config, taxid)) for taxid in taxids]
    drawn = [(taxid, species) for taxid, species in drawn if species is not None]
    if not drawn:
        return

    # every row of df is the same parasite, which is what the page selected before
    # filtering; filter_tissues in web_utils reads it the same way
    infected = infected_organs(config, df['taxid1'].unique()[0])
    if not infected:
        if title_as_subheader:
            st.subheader('Where the predicted interactions can take place')
        else:
            st.markdown('##### Where the predicted interactions can take place')
        st.caption('The figure draws none of the tissues this parasite is recorded as '
                   'infecting, so there is nothing to shade.')
        return

    # the organs the filter is set to, which are outlined on a clickable figure, and the
    # organs shaded: with nothing selected they are the same, every organ the parasite
    # infects
    filtered_organs = tissue_organs(config, selected_tissues) if selected_tissues else set()
    shown_organs = infected & (filtered_organs or infected)

    if title_as_subheader:
        st.subheader('Where the predicted interactions can take place')
    else:
        st.markdown('##### Where the predicted interactions can take place')
    st.caption('Predicted interactions whose host protein is expressed in each organ, '
               'after the confidence score and the tissue filters, and only in the organs '
               'this parasite is recorded as infecting. TISSUES annotates a host protein '
               'to every organ it is detected in, so an interaction is counted in each of '
               'the ones shown and the organs can add up to more than the network.'
               + (' Click an organ to filter the predictions on the tissues it stands for, '
                  'and click it again to clear them.' if clickable else ''))
    figures = []
    for taxid, species in drawn:
        svg, organs = load_figure(species)
        if svg is None:
            continue

        host_df = df[df['taxid2'] == str(taxid)]
        counts = count_interactions(host_df, figure_tissues)
        # An organ the annotation knows but this figure does not draw would otherwise
        # stretch the colour scale to a range nothing on the figure can reach, and an
        # organ the parasite does not infect is not somewhere the interaction can happen.
        counts = {organ: count for organ, count in counts.items()
                  if organ in organs and organ in shown_organs}
        figures.append((taxid, svg, counts, infected & organs))

    if not figures:
        return

    shared_bounds = None
    if shared_color_scale:
        shared_counts = {(taxid, organ): count
                         for taxid, _, counts, _ in figures
                         for organ, count in counts.items()}
        shared_bounds, _ = color_scale(shared_counts)

    if len(figures) > 1:
        labels = st.columns(len(figures))
        for column, (taxid, _, _, _) in zip(labels, figures):
            with column:
                st.caption(config['hosts'][int(taxid)]['label'])

    columns = st.columns(len(figures))
    for column, (taxid, svg, counts, drawn_organs) in zip(columns, figures):
        bounds, highest = (shared_bounds, max(counts.values(), default=0)
                           ) if shared_bounds is not None else color_scale(counts)

        with column:
            figure = shade_figure(svg, counts, bounds)
            figure = bottom_aligned(figure) if shared_color_scale else figure
            if compact_human and not shared_color_scale and get_species(config, taxid) == 'human':
                figure = compact_human_views(figure)
            if clickable and click_detector is not None:
                figure = link_organs(figure, drawn_organs, filtered_organs, config)
                click_detector(figure, key=f'{ORGAN_CLICK_KEY}_{taxid}')
            else:
                st.markdown(figure, unsafe_allow_html=True)
            if highest and not shared_color_scale:
                st.markdown(legend_html(bounds), unsafe_allow_html=True)
            else:
                if not highest:
                    st.caption('None of the host proteins of this network are annotated to '
                               'an organ this parasite infects.')

    if shared_color_scale and shared_bounds:
        st.markdown(legend_html(shared_bounds, compact=True), unsafe_allow_html=True)
