"""
Shades the TISSUES body figure of a host with the number of predicted interactions
reaching each of its organs.

The figures (images/tissues_<species>.svg) come from tissues.jensenlab.org and draw 21
organs, each element carrying the organ as its `title` attribute. That is a coarser
vocabulary than the 33 lifecycle tissues of the configuration, which is what decides
whether a prediction is shown at all; the organs are read from the same TISSUES download
by scripts/build_figure_tissues.py, so nothing has to be mapped between the two.
"""
import os
import xml.etree.ElementTree as ET

import pandas as pd
import streamlit as st

import utils
# the same BTO code -> organ map the annotation was built with, rather than a second copy
# of it here: the two have to name the 21 organs identically or the shading silently
# misses some
from scripts.build_figure_tissues import FIGURE_ORGANS

SVG_NS = 'http://www.w3.org/2000/svg'
ET.register_namespace('', SVG_NS)
ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')

FIGURE_DIR = 'images'

# The figures do not all label the same organ the same way: the pig says `thyroid` where
# the others say `thyroid gland`, and the annotation follows the majority.
ORGAN_ALIASES = {'thyroid': 'thyroid gland'}

# The green ramp of the TISSUES figures, as their own legend uses it. The first colour is
# for an organ no interaction reaches, so it stays the white the figures are drawn in.
NO_INTERACTIONS_COLOR = '#ffffff'
PALETTE = ['#dcf2d7', '#acdea6', '#6dc072', '#30994f', '#006d2c']

# the legend is stripped and redrawn below the figure: it is labelled with the confidence
# scores of the TISSUES website, which are not what is being shown here, and the pig
# draws those labels as outlines rather than text, so they cannot simply be rewritten
LEGEND_ID = 'Legend'

# drawing units left between the bottom of the body and the edge of the cropped figure
CROP_MARGIN = 10

# The tissues a parasite infects are the 33 fine-grained BTO terms of config['tissues'];
# the figure draws 21 coarser organs (20 for rat, which has no gall bladder). Of the 27
# terms the parasites actually use, twelve are drawn as themselves and come straight out
# of FIGURE_ORGANS. These are the five of the rest that have an organ containing them on
# the figure -- seven parasites are recorded only in `small intestine`, which the figure
# draws as `intestine`, and would otherwise shade nothing at all. (`gastrointestinal
# tract` is in config['tissues'] but no parasite uses it; it is kept here so the map
# covers the whole vocabulary.)
#
# Ten terms still have nothing standing for them:
#   `bile duct`         Clonorchis sinensis, Opisthorchis felineus, O. viverrini
#   `lymph vessel`      Brugia malayi, B. timori, Wuchereria bancrofti
#   `macrophage`        Leishmania braziliensis, L. infantum, L. major
#   `mesenteric artery` Angiostrongylus costaricensis
#   `mouth`             Leishmania braziliensis
#   `nose`              Leishmania braziliensis
#   `placenta`          Toxoplasma gondii
#   `urethra`           Trichomonas vaginalis
#   `urinary bladder`   Schistosoma haematobium
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
}


@st.cache_data(show_spinner=False)
def load_figure_tissues(data_dir):
    '''
    Host proteins annotated with the organs of the body figure, written by
    scripts/build_figure_tissues.py. The file is not part of the older snapshot data
    directories, so a missing one only leaves the figure out.

    :param str data_dir: directory holding figure_tissues.parquet
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
    both the download and the figure, images/tissues_<species>.svg.

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
    expressed in one of those tissues, but TISSUES then annotates that protein to about
    three organs, most of which the parasite never reaches: a Loa loa protein selected for
    being expressed in skin also comes annotated to the nervous system, and the figure drew
    the brain as the darkest organ on the page. Restricting the shading to the organs the
    parasite actually infects is what keeps the figure about the parasite rather than about
    how broadly its targets happen to be expressed.

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


def count_interactions(df, figure_tissues):
    '''
    Counts the predicted interactions reaching each organ. An interaction is counted once
    per organ its host protein is annotated to, so the counts do not add up to the size
    of the network -- a host protein is annotated to about three organs.

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
    # the busiest organ is the darkest green whatever the size of the network
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


def legend_html(bounds):
    '''The colour bins as a row of swatches, labelled with the counts they stand for.'''
    swatches = []
    previous = 0
    for upper, color in bounds:
        label = str(upper) if upper == previous + 1 else f'{previous + 1}-{upper}'
        swatches.append(
            f'<div style="text-align: center; flex: 1;">'
            f'<div style="background: {color}; border: 1px solid #939598; height: 14px;"></div>'
            f'<div style="font-size: 0.7rem; color: #555;">{label}</div></div>')
        previous = upper

    return ('<div style="display: flex; gap: 2px; margin-top: 0.5rem;">'
            + ''.join(swatches) + '</div>')


def show_body_figure(config, data_dir, df, taxids):
    '''
    Draws the body figure of each selected host, its organs shaded by the number of
    predicted interactions reaching them. A host group covering more than one species
    (Rodent is rat and mouse) gets one figure each, since the two are annotated
    separately and TISSUES draws a different body for each.

    :param dict config: parsed configuration
    :param str data_dir: directory holding figure_tissues.parquet
    :param df: predictions dataframe, already filtered to what the network shows
    :param taxids: taxids of the selected host
    '''
    figure_tissues = load_figure_tissues(data_dir)
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
        st.markdown('##### Where the predicted interactions can take place')
        st.caption('The figure draws none of the tissues this parasite is recorded as '
                   'infecting, so there is nothing to shade.')
        return

    st.markdown('##### Where the predicted interactions can take place')
    columns = st.columns(len(drawn))
    for column, (taxid, species) in zip(columns, drawn):
        svg, organs = load_figure(species)
        if svg is None:
            continue

        host_df = df[df['taxid2'] == str(taxid)]
        counts = count_interactions(host_df, figure_tissues)
        # an organ the annotation knows but this figure does not draw would otherwise
        # stretch the colour scale to a range nothing on the figure can reach, and an
        # organ the parasite does not infect is not somewhere the interaction can happen
        counts = {organ: count for organ, count in counts.items()
                  if organ in organs and organ in infected}
        bounds, highest = color_scale(counts)

        with column:
            if len(drawn) > 1:
                st.caption(config['hosts'][int(taxid)]['label'])
            st.markdown(shade_figure(svg, counts, bounds), unsafe_allow_html=True)
            if highest:
                st.markdown(legend_html(bounds), unsafe_allow_html=True)
            else:
                st.caption('None of the host proteins of this network are annotated to '
                           'an organ this parasite infects.')

    st.caption('Predicted interactions whose host protein is expressed in each organ, '
               'after the confidence score and the tissue filters, and only in the organs '
               'this parasite is recorded as infecting. TISSUES annotates a host protein '
               'to about three organs, so an interaction is counted in each of the ones '
               'shown and the organs can add up to more than the network.')
