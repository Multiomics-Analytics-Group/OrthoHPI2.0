"""
Builds the figures of the OrthoHPI 2.0 poster, at poster scale rather than screen scale.

The app draws these for a browser, and neither the type nor the circos survives being
printed at A0: the labels come out around 4pt on an 841mm sheet, and the chords of the
circos only take their colour when a parasite is hovered, so on paper the circle is a grey
wash. The figures are rebuilt here from the app's own data functions -- so they still say
what the app says -- with the type sized for reading at a distance and the circos coloured
statically.

The hand-drawn figures are written in millimetres: their viewBox is their final size on
the poster, so a font-size in the file is a font-size on the printed sheet and nothing has
to be back-calculated through a scale factor. The two figures that come from elsewhere
(the Plotly bar chart and the body figure) carry their own coordinates and are scaled into
place by the assembler, so their type is set through px_for_pt instead.

Writes to poster/figures/. Nothing in app/, config.yml or the pipeline is touched.

    .venv/bin/python scripts/build_poster_figures.py
"""
import importlib.util
import logging
import math
import os
import subprocess
import sys
import types
import xml.etree.ElementTree as ET

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'app'))

# the app pages are written to run inside streamlit; imported bare they warn about the
# missing script run context on every cached call, which is noise here and nothing else
logging.getLogger('streamlit').setLevel(logging.ERROR)

import pandas as pd  # noqa: E402
from PIL import Image  # noqa: E402

import utils  # noqa: E402
import web_utils  # noqa: E402
import body_figure  # noqa: E402
import structure_visualizer as strv  # noqa: E402

DATA_DIR = os.path.join(REPO, 'data')
CONFIG_FILE = os.path.join(REPO, 'config.yml')
OUT_DIR = os.path.join(REPO, 'poster', 'figures')
ASSET_DIR = os.path.join(REPO, 'poster', 'assets')

MM_PER_PT = 25.4 / 72

# The final size of each figure on the poster, in millimetres. The hand-drawn figures are
# written straight into these coordinates; the assembler places every figure in a box of
# exactly this size, so changing a number here moves the poster with it.
SIZES = {
    'interactions': (780, 130),
    'shared': (780, 155),
    'circos': (380, 290),
    'network': (395, 290),
    # the TISSUES drawing is a front and a side view side by side, so it comes out wider
    # than tall; the figure is cropped to its content and this is the box that shape is
    # fitted into, not a shape imposed on it
    'body': (170, 105),
}

# The example interactome the poster zooms into, and the one interaction of it the last
# panel opens up. The pair is predicted at 0.944 and sits inside the drawn network, so the
# last two steps of the funnel are the same interaction seen at two scales.
EXAMPLE_PARASITE = 'Plasmodium falciparum'
EXAMPLE_HOST = 'Homo sapiens'
EXAMPLE_PAIR = [
    # (panel, species, protein, UniProt, description)
    ('parasite', 'Plasmodium falciparum', 'PF3D7_1104100', 'Q8IIW1', 'Synaptobrevin, putative'),
    ('host', 'Homo sapiens', 'YKT6', 'O15498', 'Synaptobrevin homolog YKT6'),
]

# Confidence floor for the example network. The predictions are not restricted to the
# tissues the parasite infects here -- P. falciparum is drawn as the whole interactome --
# and unfiltered it is 297 interactions over 222 proteins, which no poster panel can label.
# 0.6 leaves 42 proteins, about what the panel carried before, and keeps the pair the last
# panel opens up (0.944) inside the picture.
NETWORK_SCORE = 0.6
# components smaller than this are left out of the drawn network -- see build_network
NETWORK_MIN_COMPONENT = 3

# Host proteins shown in the dot matrix, most-shared first. The figure is bounded by its
# rows: a row needs about 9mm to carry an 18pt protein label, and the band can spare 155mm
# once the x axis has its 31 rotated species names and the legend has its strip.
SHARED_TOP = 12

# chords are kept only for the parasite pairs sharing at least this many host proteins.
# Colouring all 377 chords of the human view at once is what the app avoids by colouring
# on hover -- see generate_circos_plot in app/pages/2_Compare_Parasites.py -- because the
# circle becomes a wash no scale can be read off. The threshold is reported when the
# figure is built so it can be retuned against the number of chords it leaves.
CHORD_MIN_SHARED = 10

# the ramp the app already carries the shared-protein count on, so the poster and the app
# colour the same quantity the same way
CHORD_CMAP = ['#a6bddb', '#74a9cf', '#3690c0', '#0570b0', '#045a8d', '#034368', '#023858']

INK = '#1d2126'
MUTED = '#6b7280'
RULE = '#d4d8dd'
UNKNOWN_COLOR = '#999999'


def px_for_pt(points, export_px_width, target_mm_width):
    """
    The font size, in the pixels of an exported figure, that lands at `points` once that
    figure has been scaled into a box `target_mm_width` wide on the printed poster.

    A figure that carries its own coordinate system -- anything Plotly writes -- is scaled
    bodily into its box, so a size set in it is only meaningful relative to the width of
    the export. This converts the size that is actually wanted on paper into that.
    """
    return points * MM_PER_PT * export_px_width / target_mm_width


def load_app_page(relative_path, name):
    """
    Imports one of the streamlit pages for its data functions.

    The pages run their body on import -- headings, captions, charts -- which is harmless
    outside a running app, every st.* call becoming a no-op that warns. Only
    streamlit_bokeh refuses to import at all, registering a component that needs a real
    app behind it, so it is stubbed before the page reaches it. Nothing the poster uses
    goes near bokeh: the circos is redrawn here from the dataframes.
    """
    if 'streamlit_bokeh' not in sys.modules:
        stub = types.ModuleType('streamlit_bokeh')
        stub.streamlit_bokeh = lambda *args, **kwargs: None
        sys.modules['streamlit_bokeh'] = stub

    spec = importlib.util.spec_from_file_location(name, os.path.join(REPO, relative_path))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)

    return module


def svg_header(width_mm, height_mm):
    """Opens an SVG whose user unit is one millimetre of the printed poster."""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{width_mm}mm" height="{height_mm}mm" '
            f'viewBox="0 0 {width_mm} {height_mm}">\n')


def write_svg(name, body):
    path = os.path.join(OUT_DIR, name)
    with open(path, 'w') as handle:
        handle.write(body)
    print(f'  wrote {name} ({os.path.getsize(path) / 1024:.0f} KB)')

    return path


def escape(text):
    return (str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;'))


def crop_to_drawing(path):
    """
    Shrinks an SVG's viewBox onto what it actually draws, through Inkscape.

    Only worth doing for a figure that arrives with a viewBox larger than its content;
    a bounding box over arbitrary path data is not something to work out by hand, and
    Inkscape is already what the poster is edited in. If it is not installed the figure
    is simply left as it was, which costs some white space and nothing else.
    """
    inkscape = '/Applications/Inkscape.app/Contents/MacOS/inkscape'
    if not os.path.exists(inkscape):
        print(f'    (Inkscape not found, {os.path.basename(path)} left uncropped)')
        return

    result = subprocess.run(
        [inkscape, '--export-area-drawing', '--export-plain-svg',
         f'--export-filename={path}', path],
        capture_output=True, text=True)
    if result.returncode != 0:
        print(f'    (cropping {os.path.basename(path)} failed: {result.stderr.strip()[:120]})')


# ----------------------------------------------------------------------------------
# 1. how many interactions each parasite has -- the scale of the resource
# ----------------------------------------------------------------------------------

def build_interactions_per_parasite(config):
    """
    The app's own bar chart, restyled for print: the type set through px_for_pt so it
    lands at a readable size once scaled onto the poster, and the legend moved inside the
    plotting area, where on a poster there is no page width to spend on it.
    """
    home = load_app_page('app/OrthoHPI_Home.py', 'poster_home')

    overview = home.get_overview_predictions(DATA_DIR, config)
    figure = home.generate_interactions_per_parasite(overview, config.get('parasite_groups', {}))

    width_mm, height_mm = SIZES['interactions']
    export_px = 2000
    tick_px = px_for_pt(15, export_px, width_mm)
    axis_px = px_for_pt(18, export_px, width_mm)
    legend_px = px_for_pt(17, export_px, width_mm)

    # the margins are what the figure spends on its labels, so they belong in the
    # millimetres those labels are read at, not in pixels of an export whose size is
    # arbitrary: as fixed pixels they swallowed the whole plotting area the moment the
    # figure was made shorter
    per_mm = export_px / width_mm
    figure.update_layout(
        width=export_px, height=round(export_px * height_mm / width_mm),
        # the top margin carries two rows -- the legend, and the host name over each
        # column that make_subplots writes as an annotation -- and the bottom one the
        # forty species names, set on their side
        margin=dict(l=round(26 * per_mm), r=round(4 * per_mm),
                    t=round(26 * per_mm), b=round(50 * per_mm)),
        font=dict(family='Arial', color=INK),
        legend=dict(orientation='h', yanchor='bottom', y=1.10, x=0,
                    font=dict(size=legend_px), title_text=''),
        plot_bgcolor='white', paper_bgcolor='white')
    figure.update_xaxes(tickfont=dict(size=tick_px), tickangle=-60,
                        showline=True, linecolor=RULE, linewidth=2)
    figure.update_yaxes(tickfont=dict(size=tick_px), title_font=dict(size=axis_px),
                        gridcolor='#eef0f2', showline=False)
    # the host name over each column, which make_subplots writes as annotations
    for annotation in figure.layout.annotations:
        annotation.font.size = axis_px
        annotation.font.color = INK

    path = os.path.join(OUT_DIR, 'fig1-interactions-per-parasite.svg')
    figure.write_image(path, format='svg')
    print(f'  wrote fig1-interactions-per-parasite.svg '
          f'({os.path.getsize(path) / 1024:.0f} KB, {overview["taxid1_label"].nunique()} parasites)')

    return overview


# ----------------------------------------------------------------------------------
# 1b. what the parasites share -- the host proteins many of them reach
# ----------------------------------------------------------------------------------

def build_shared_proteins(config, compare):
    """
    The host proteins the most parasites are predicted to reach, as a dot per pair.

    The circos says which parasites overlap; nothing else on the poster says what they
    overlap on, which is the question a reader of the circos asks next. Each dot is sized
    by how many of that parasite's proteins reach that host protein, so a dot is not merely
    a yes.
    """
    predictions = compare.get_tissue_expressed_predictions(
        DATA_DIR, config,
        tuple(web_utils.get_host_groups(
            config, web_utils.load_predictions(DATA_DIR))[EXAMPLE_HOST]))

    groups = {p['label']: p.get('group', 'Unclassified') for p in config['parasites'].values()}
    palette = config.get('parasite_groups', {})
    annotations = web_utils.load_protein_annotations(DATA_DIR)

    result = compare.get_top_shared_proteins(
        predictions, groups, {g: i for i, g in enumerate(palette)}, annotations,
        top=SHARED_TOP)
    if result is None:
        raise SystemExit('no host protein is shared by more than one parasite')
    dots, proteins, parasites, most = result

    figure = compare.generate_shared_protein_dots(dots, proteins, parasites, most, palette)

    width_mm, height_mm = SIZES['shared']
    export_px = 2000
    per_mm = export_px / width_mm
    tick_px = px_for_pt(18, export_px, width_mm)
    legend_px = px_for_pt(18, export_px, width_mm)

    figure.update_layout(
        width=export_px, height=round(export_px * height_mm / width_mm),
        # the left margin is left to plotly: the row labels are a gene symbol and a
        # description, and how wide that runs is not something to guess at here
        margin=dict(l=0, r=round(4 * per_mm), t=round(14 * per_mm), b=round(46 * per_mm)),
        font=dict(family='Arial', color=INK),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, x=0,
                    font=dict(size=legend_px), title_text=''),
        yaxis_title=None, plot_bgcolor='white', paper_bgcolor='white')
    figure.update_xaxes(tickfont=dict(size=tick_px), tickangle=-60,
                        showline=True, linecolor=RULE, linewidth=2, gridcolor='#eef0f2')
    # dtick=1 forces a label on every row. Left to itself plotly thins the ticks when it
    # decides they are too close for the height, and drops every second protein -- half the
    # point of the figure is that the rows are named
    figure.update_yaxes(tickfont=dict(size=tick_px), gridcolor='#eef0f2', automargin=True,
                        dtick=1, tick0=0)
    # the dots have to grow with the export or they are specks on a 2000px canvas
    figure.update_traces(marker=dict(sizemin=round(2.2 * per_mm)))
    figure.update_layout(**{'yaxis_ticksuffix': '  '})

    path = os.path.join(OUT_DIR, 'fig1b-shared-proteins.svg')
    figure.write_image(path, format='svg')
    print(f'  wrote fig1b-shared-proteins.svg ({os.path.getsize(path) / 1024:.0f} KB, '
          f'{len(proteins)} proteins x {len(parasites)} parasites, '
          f'most-shared reaches {most} parasites)')


# ----------------------------------------------------------------------------------
# 2. which parasites reach the same host proteins -- the circos
# ----------------------------------------------------------------------------------

def chord_path(a1, a2, radius):
    """
    One chord, as a quadratic bezier between two points on the circle with the centre as
    its control point. The curve is then pulled towards the middle by an amount that
    follows how far apart its ends are: a pair of neighbours gets a shallow arc hugging
    the rim, a pair facing each other gets something close to a diameter, which is what
    keeps the middle of a circle this busy from filling in solid.
    """
    x1, y1 = radius * math.cos(a1), radius * math.sin(a1)
    x2, y2 = radius * math.cos(a2), radius * math.sin(a2)

    separation = abs(((a2 - a1) + math.pi) % (2 * math.pi) - math.pi) / math.pi
    pull = 1 - separation  # 0 for opposite ends of the circle, ~1 for neighbours
    mid = ((a1 + a2) / 2) if abs(a2 - a1) <= math.pi else ((a1 + a2) / 2 + math.pi)
    cx, cy = pull * radius * math.cos(mid), pull * radius * math.sin(mid)

    return f'M{x1:.2f},{y1:.2f} Q{cx:.2f},{cy:.2f} {x2:.2f},{y2:.2f}'


def chord_color(shared, lowest, highest):
    """The shared-protein count on the app's ramp, logarithmically: half the pairs share
    fewer than five proteins and one shares thirty-odd, so a linear scale puts almost
    every chord in the first colour."""
    if highest <= lowest:
        return CHORD_CMAP[-1]
    position = (math.log(shared) - math.log(lowest)) / (math.log(highest) - math.log(lowest))

    return CHORD_CMAP[min(int(position * len(CHORD_CMAP)), len(CHORD_CMAP) - 1)]


def build_circos(config, compare):
    """
    The parasites of the human view around a circle, joined where they are predicted to
    reach the same host proteins.

    Only the pairs sharing at least CHORD_MIN_SHARED proteins are drawn, and those carry
    the count as colour. The app instead draws every pair grey and colours one parasite's
    chords on hover; printed, that is a grey wash, and drawing all of them coloured is the
    wash the app's own comment warns about. Cutting the thin end keeps the figure about
    the strong overlaps, which is the thing worth saying at poster distance.
    """
    predictions = compare.get_tissue_expressed_predictions(
        DATA_DIR, config, tuple(web_utils.get_host_groups(config, web_utils.load_predictions(DATA_DIR))[EXAMPLE_HOST]))

    groups = {p['label']: p.get('group', 'Unclassified') for p in config['parasites'].values()}
    palette = config.get('parasite_groups', {})
    links, nodes = compare.get_common_interactors(predictions, groups,
                                                  {g: i for i, g in enumerate(palette)})

    kept = links[links['shared'] >= CHORD_MIN_SHARED]
    print(f'  circos: {len(nodes)} parasites, {len(links)} sharing pairs, '
          f'{len(kept)} kept at >= {CHORD_MIN_SHARED} shared proteins')

    width, height = SIZES['circos']
    # the circle is sized from the box rather than fixed, because what constrains it is
    # not the circle but the species names standing off it: a label runs outward along its
    # own radius, so the drawing is a good deal wider than the circle in every direction
    LEGEND_SPACE = 34      # the clade key and the chord ramp along the bottom
    LABEL_SIZE = 6.4
    LABEL_ALLOWANCE = 44   # the longest of the abbreviated species names, set at that size
    cx, cy = width / 2, (height - LEGEND_SPACE) / 2
    radius = min(cy, width / 2) - LABEL_ALLOWANCE - 16
    arc_radius = radius + 6      # the rim carrying the group colour
    label_radius = radius + 12

    lowest, highest = kept['shared'].min(), kept['shared'].max()
    count = len(nodes)
    step = 2 * math.pi / count
    # the circle starts at the top and runs clockwise, which is the order the labels are
    # read in and the order the app lists the parasites in
    angles = {i: -math.pi / 2 + i * step for i in range(count)}

    out = [svg_header(width, height)]
    out.append(f'<g transform="translate({cx},{cy})" font-family="Arial">\n')

    # chords first, so the rim and the labels sit over them
    out.append('<g fill="none" stroke-linecap="round">\n')
    for _, row in kept.sort_values('shared').iterrows():
        color = chord_color(row['shared'], lowest, highest)
        # the strongest overlaps are drawn a little heavier as well as darker: colour
        # alone is hard to rank at three metres
        weight = 0.35 + 0.5 * (row['shared'] - lowest) / max(highest - lowest, 1)
        out.append(f'<path d="{chord_path(angles[row["source"]], angles[row["target"]], radius)}" '
                   f'stroke="{color}" stroke-width="{weight:.2f}" opacity="0.75"/>\n')
    out.append('</g>\n')

    # the rim, one arc per parasite, coloured by taxonomic group
    gap = step * 0.16
    for _, node in nodes.iterrows():
        angle = angles[node['index']]
        start, end = angle - step / 2 + gap / 2, angle + step / 2 - gap / 2
        x1, y1 = arc_radius * math.cos(start), arc_radius * math.sin(start)
        x2, y2 = arc_radius * math.cos(end), arc_radius * math.sin(end)
        color = palette.get(node['group'], UNKNOWN_COLOR)
        out.append(f'<path d="M{x1:.2f},{y1:.2f} A{arc_radius},{arc_radius} 0 0 1 {x2:.2f},{y2:.2f}" '
                   f'fill="none" stroke="{color}" stroke-width="4.5"/>\n')

        # the label runs outward along its own radius, flipped on the left hand side so it
        # is never upside down
        degrees = math.degrees(angle)
        flip = 90 < (degrees % 360) < 270
        lx, ly = label_radius * math.cos(angle), label_radius * math.sin(angle)
        rotate = degrees + 180 if flip else degrees
        anchor = 'end' if flip else 'start'
        out.append(f'<text x="{lx:.2f}" y="{ly:.2f}" font-size="{LABEL_SIZE}" fill="{INK}" '
                   f'text-anchor="{anchor}" dominant-baseline="middle" '
                   f'transform="rotate({rotate:.2f} {lx:.2f} {ly:.2f})">'
                   f'<tspan font-style="italic">{escape(node["name"])}</tspan></text>\n')

    out.append('</g>\n')
    out.append(circos_scales(width, height, palette, set(nodes['group']), lowest, highest))
    out.append('</svg>\n')

    write_svg('fig2-circos.svg', ''.join(out))


def circos_scales(width, height, palette, shown, lowest, highest):
    """The two keys the circle needs: which colour is which clade, and what the darkness
    of a chord means. v1 carried neither."""
    out = ['<g font-family="Arial">\n']

    # taxonomic groups, along the bottom in two rows
    entries = [(g, c) for g, c in palette.items() if g in shown]
    per_row = math.ceil(len(entries) / 2)
    for i, (group, color) in enumerate(entries):
        col, row = i % per_row, i // per_row
        x = 6 + col * (width - 90) / per_row
        y = height - 20 + row * 8
        out.append(f'<rect x="{x}" y="{y - 3.4}" width="5" height="5" rx="1" fill="{color}"/>\n')
        out.append(f'<text x="{x + 7.5}" y="{y}" font-size="6.4" fill="{INK}" '
                   f'dominant-baseline="middle">{escape(group)}</text>\n')

    # the chord ramp, upright in the bottom right corner
    bar_x, bar_y, bar_w, bar_h = width - 58, height - 24, 46, 5
    step = bar_w / len(CHORD_CMAP)
    for i, color in enumerate(CHORD_CMAP):
        out.append(f'<rect x="{bar_x + i * step:.2f}" y="{bar_y}" width="{step + 0.3:.2f}" '
                   f'height="{bar_h}" fill="{color}"/>\n')
    out.append(f'<text x="{bar_x}" y="{bar_y - 3}" font-size="6" fill="{MUTED}">'
               f'shared host proteins</text>\n')
    out.append(f'<text x="{bar_x}" y="{bar_y + bar_h + 6.4}" font-size="6" fill="{MUTED}">'
               f'{lowest}</text>\n')
    out.append(f'<text x="{bar_x + bar_w}" y="{bar_y + bar_h + 6.4}" font-size="6" '
               f'fill="{MUTED}" text-anchor="end">{highest}</text>\n')
    out.append('</g>\n')

    return ''.join(out)


# ----------------------------------------------------------------------------------
# 3. one interactome: where it acts, and what it touches
# ----------------------------------------------------------------------------------

def example_predictions(config):
    """
    The predictions the panels of the example interactome are drawn from: one parasite in
    one host, and every interaction predicted for it.

    The app's own page restricts these to the tissues the parasite infects, and so did an
    earlier version of the poster. The network panel shows the whole interactome instead:
    the tissue restriction is a statement about where an interaction could happen, which
    the body figure beside it already makes, and applying it twice left the network looking
    like the whole prediction when it was a filtered part of it.
    """
    host_taxids = web_utils.get_host_groups(
        config, web_utils.load_predictions(DATA_DIR))[EXAMPLE_HOST]
    df = web_utils.get_host_predictions(DATA_DIR, tuple(host_taxids))

    return df[df['taxid1_label'] == EXAMPLE_PARASITE]


def build_body_figure(config, df):
    """The TISSUES body figure with its organs shaded by how many predicted interactions
    reach them. Already an SVG and already assembled by the app, so it is only re-emitted
    at a fixed size with the swatch key drawn beside it."""
    figure_tissues = body_figure.load_figure_tissues(DATA_DIR)
    svg, organs = body_figure.load_figure('human')

    # the organs this parasite is recorded as infecting, the same restriction the app
    # applies. Without it the figure shades wherever the targeted host proteins happen to
    # be expressed, which for a protein annotated to about three organs is most of the body
    infected = body_figure.infected_organs(config, df['taxid1'].unique()[0])

    counts = body_figure.count_interactions(df, figure_tissues)
    counts = {organ: n for organ, n in counts.items()
              if organ in organs and organ in infected}
    bounds, highest = body_figure.color_scale(counts)
    shaded = body_figure.shade_figure(svg, counts, bounds)

    root = ET.fromstring(shaded)
    root.attrib.pop('style', None)
    # width and height have to agree with the viewBox or the drawing is stretched, and a
    # stretched drawing is also what would then be measured when it is cropped. The
    # assembler scales the figure into its box afterwards, so these are only ever a shape.
    x, y, view_width, view_height = (
        float(v) for v in root.get('viewBox').replace(',', ' ').split())
    root.set('width', str(view_width))
    root.set('height', str(view_height))

    path = write_svg('fig3a-body.svg', ET.tostring(root, encoding='unicode'))
    # the source drawing sits in a viewBox much taller than itself -- the app lets the
    # browser absorb that, but inlined on the poster it would be a band of empty space
    # above and below the body
    crop_to_drawing(path)

    cropped = ET.parse(path).getroot()
    _, _, final_width, final_height = (
        float(v) for v in cropped.get('viewBox').replace(',', ' ').split())
    print(f'    organs shaded: {len(counts)}, busiest {highest} interactions, '
          f'cropped aspect {final_width / final_height:.2f}')

    return bounds


def pack_components(graph, components, width, height):
    """
    Lays each connected component out on its own and packs them into the panel.

    The predictions of one parasite do not make one connected network: Loa loa comes out
    as a couple of hubs and a tail of isolated pairs. Handed the whole graph, spring_layout
    has nothing to push the pieces apart with and drops them on top of each other, so each
    piece is laid out separately and given a cell of its own, with the area of a cell
    following the number of proteins in it. The cells are shelf-packed -- placed left to
    right, wrapping to a new row when the next one will not fit -- which for a handful of
    boxes of similar size wastes little enough space and keeps the big hubs at the top,
    where they are read first.

    :return: {node: (x, y)} in millimetres within the panel
    """
    import networkx as nx

    CELL_UNIT = 13      # cell side per protein**0.62, before the packing is scaled to fit
    CELL_FLOOR = 24     # a component still needs room for its labels
    PAD = 4             # between cells, in the same units

    # area grows a little slower than the protein count so that a hub, whose labels sit
    # around its rim rather than through the middle, is not given the whole panel
    boxes = [(c, max(len(c) ** 0.62 * CELL_UNIT, CELL_FLOOR)) for c in components]

    def shelf_pack(shelf_width):
        placed, x, shelf_y, shelf_height = [], 0.0, 0.0, 0.0
        for component, side in boxes:
            if x > 0 and x + side > shelf_width:
                x, shelf_y, shelf_height = 0.0, shelf_y + shelf_height + PAD, 0.0
            placed.append((component, side, x, shelf_y))
            x += side + PAD
            shelf_height = max(shelf_height, side)
        # the right edge of the widest row, not the left edge of its last cell
        used = max((cell_x + side for _, side, cell_x, _ in placed), default=shelf_width)

        return placed, used, shelf_y + shelf_height

    # A shelf width has to be chosen before anything can be packed, and the one that fills
    # the panel best depends on the mix of component sizes. Rather than guess a multiplier,
    # every width that could start a new row is tried and the one the panel can be scaled
    # up the most for wins.
    margin = 12
    candidates = sorted({sum(s for _, s in boxes[:i]) + PAD * i
                         for i in range(1, len(boxes) + 1)} | {max(s for _, s in boxes)})
    best = None
    for candidate in candidates:
        placed, packed_width, packed_height = shelf_pack(candidate)
        scale = min((width - 2 * margin) / packed_width, (height - 2 * margin) / packed_height)
        if best is None or scale > best[0]:
            best = (scale, placed, packed_width, packed_height)

    scale, placed, packed_width, packed_height = best
    # centre the packed block in the panel
    offset_x = (width - packed_width * scale) / 2
    offset_y = (height - packed_height * scale) / 2

    pos = {}
    for component, side, cell_x, cell_y in placed:
        subgraph = graph.subgraph(component)
        # a fixed seed so the figure is the same every time it is built
        # a generous k: most of these components are stars, and at the default the leaves
        # of an eighteen-protein hub land in a ring too tight to label
        local = nx.spring_layout(subgraph, seed=11, iterations=500,
                                 k=3.4 / math.sqrt(max(len(component), 2)))
        xs = [p[0] for p in local.values()]
        ys = [p[1] for p in local.values()]
        span_x = (max(xs) - min(xs)) or 1.0
        span_y = (max(ys) - min(ys)) or 1.0
        # inset so the nodes on the edge of a component keep their labels inside its cell
        inset = side * 0.16
        usable = side - 2 * inset
        for node, (nx_, ny_) in local.items():
            fx = (nx_ - min(xs)) / span_x if span_x else 0.5
            fy = (ny_ - min(ys)) / span_y if span_y else 0.5
            pos[node] = (offset_x + (cell_x + inset + fx * usable) * scale,
                         offset_y + (cell_y + inset + fy * usable) * scale)

    return pos


def separate_labels(labels, font_size, width, height, obstacles=(), iterations=200):
    """
    Nudges overlapping protein labels apart.

    Placing a label on the far side of its node from the middle of its component sorts out
    most of them, but a hub with a dozen leaves still stacks several in the same place, and
    an unreadable name is worse than a slightly displaced one. Overlapping pairs are pushed
    apart along whichever axis they overlap least on -- usually the vertical, which keeps a
    label near the node it belongs to -- and each label is pulled back towards its own node
    afterwards so the whole set does not drift outwards.

    :param labels: [{'x','y','text','anchor'}], moved in place; 'anchor' is its node
    :param float font_size: in the same millimetres as the coordinates
    """
    for label in labels:
        # Arial's average advance is a little over half its point size; only the ratio
        # matters here, since this is deciding overlaps rather than typesetting
        label['w'] = 0.55 * font_size * len(label['text'])
        label['h'] = font_size * 1.1

    for _ in range(iterations):
        moved = False

        # the markers do not move, so a label overlapping one is pushed clear of it. Two
        # nodes drawn almost on top of each other -- ITGB1 and ITGB5 in the Loa loa hub --
        # otherwise get labels that separate from each other and land on the markers
        for label in labels:
            for ox, oy, oradius in obstacles:
                dy = abs(label['y'] - oy) - (label['h'] / 2 + oradius)
                dx = abs(label['x'] - ox) - (label['w'] / 2 + oradius)
                if dx >= 0 or dy >= 0:
                    continue
                moved = True
                label['y'] += (-dy + 0.1) * (1 if label['y'] > oy else -1)

        for i, a in enumerate(labels):
            for b in labels[i + 1:]:
                dx = abs(a['x'] - b['x']) - (a['w'] + b['w']) / 2
                dy = abs(a['y'] - b['y']) - (a['h'] + b['h']) / 2
                if dx >= 0 or dy >= 0:
                    continue
                moved = True
                if -dy <= -dx:
                    shift = (-dy) / 2 + 0.05
                    sign = 1 if a['y'] > b['y'] else -1
                    a['y'] += shift * sign
                    b['y'] -= shift * sign
                else:
                    shift = (-dx) / 2 + 0.05
                    sign = 1 if a['x'] > b['x'] else -1
                    a['x'] += shift * sign
                    b['x'] -= shift * sign

        # back towards the node, so a label separated from its neighbours does not end up
        # closer to somebody else's
        for label in labels:
            ax, ay = label['anchor']
            label['x'] += (ax - label['x']) * 0.08
            label['y'] += (ay - label['y']) * 0.04
            label['x'] = min(max(label['x'], label['w'] / 2 + 1), width - label['w'] / 2 - 1)
            label['y'] = min(max(label['y'], font_size), height - 12)

        if not moved:
            break


def build_network(config, df):
    """
    The example interactome as a drawn network.

    pyvis renders through vis.js into a canvas, which cannot be exported as vector at all,
    so the layout is recomputed with networkx and the result written as SVG. The parasite
    and host proteins keep the shapes and colours the app gives them -- a diamond against a
    circle -- so a reader who has seen the app recognises the picture.
    """
    import networkx as nx

    edges = df[df['weight'] >= NETWORK_SCORE].drop_duplicates(subset=['source', 'target'])
    graph = nx.Graph()
    for row in edges.itertuples():
        graph.add_node(row.source, name=row.source_name, kind='parasite')
        graph.add_node(row.target, name=row.target_name, kind='host')
        graph.add_edge(row.source, row.target, weight=float(row.weight))

    components = sorted(nx.connected_components(graph), key=len, reverse=True)
    # the isolated pairs are left out. A parasite protein joined to one host protein and
    # nothing else says only that the pair was predicted, which the table already says,
    # and there are eight of them: drawn, they take half the panel from the hubs that
    # carry the shape of the interactome
    components = [c for c in components if len(c) >= NETWORK_MIN_COMPONENT]
    drawn = {n for c in components for n in c}
    graph = graph.subgraph(drawn).copy()
    print(f'  network: {graph.number_of_nodes()} proteins, {graph.number_of_edges()} '
          f'interactions at score >= {NETWORK_SCORE}, '
          f'{len(components)} components {[len(c) for c in components]}')

    width, height = SIZES['network']
    pos = pack_components(graph, components, width, height)
    degrees = dict(graph.degree())
    busiest = max(degrees.values(), default=1)

    # where each component sits, so a label can be pushed away from the middle of its own
    # cluster instead of landing on a neighbour
    centroids = {}
    for component in components:
        cx = sum(pos[n][0] for n in component) / len(component)
        cy = sum(pos[n][1] for n in component) / len(component)
        for node in component:
            centroids[node] = (cx, cy)

    out = [svg_header(width, height), '<g font-family="Arial">\n']

    out.append('<g stroke="#9aa3ad" fill="none">\n')
    for a, b in graph.edges():
        (x1, y1), (x2, y2) = pos[a], pos[b]
        out.append(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
                   f'stroke-width="0.55"/>\n')
    out.append('</g>\n')

    LABEL_SIZE = 6.0
    labels = []
    markers = []
    for node, data in graph.nodes(data=True):
        x, y = pos[node]
        # a hub is drawn larger, the same thing the app does with its node sizes
        size = 2.4 + 2.2 * (degrees[node] / busiest) ** 0.6
        markers.append((x, y, size + 0.6))
        if data['kind'] == 'parasite':
            out.append(f'<rect x="{x - size:.2f}" y="{y - size:.2f}" width="{2 * size:.2f}" '
                       f'height="{2 * size:.2f}" fill="#D55E00" stroke="white" '
                       f'stroke-width="0.5" transform="rotate(45 {x:.2f} {y:.2f})"/>\n')
        else:
            out.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{size:.2f}" fill="#3f4650" '
                       f'stroke="white" stroke-width="0.5"/>\n')

        # the label starts on the far side of the node from the middle of its component,
        # so the labels of a cluster fan outwards rather than piling over its centre
        above = y <= centroids[node][1]
        ly = y - size - 1.8 if above else y + size + 4.4
        labels.append({'x': x, 'y': ly, 'text': str(data['name']), 'anchor': (x, ly)})

    separate_labels(labels, LABEL_SIZE, width, height, obstacles=markers)
    for label in labels:
        out.append(f'<text x="{label["x"]:.2f}" y="{label["y"]:.2f}" font-size="{LABEL_SIZE}" '
                   f'fill="{INK}" text-anchor="middle">{escape(label["text"])}</text>\n')

    out.append(f'<g font-size="6.4" fill="{MUTED}">\n')
    out.append(f'<rect x="4" y="{height - 9.6}" width="4.4" height="4.4" fill="#D55E00" '
               f'transform="rotate(45 6.2 {height - 7.4})"/>\n')
    out.append(f'<text x="12" y="{height - 6}" dominant-baseline="middle">parasite protein</text>\n')
    out.append(f'<circle cx="64" cy="{height - 7.4}" r="2.4" fill="#3f4650"/>\n')
    out.append(f'<text x="70" y="{height - 6}" dominant-baseline="middle">host protein</text>\n')
    out.append('</g>\n')

    out.append('</g>\n</svg>\n')
    write_svg('fig3b-network.svg', ''.join(out))


# ----------------------------------------------------------------------------------
# 4. the interacting pair, down to structure
# ----------------------------------------------------------------------------------

# The AlphaFold models are screenshots of the app's own 3Dmol viewer, saved by hand into
# images/ and embedded by the assembler. A backbone trace was drawn here as vector for a
# while -- there is no molecular renderer on this machine and 3Dmol only draws into a WebGL
# canvas -- but a trace has no ribbons on its sheets and no cylinders on its helices, and a
# real cartoon at a slightly soft resolution reads better on a poster than a sharp
# approximation of one. The models are still fetched, so the poster and the app show the
# same file and the download stays cached.

def build_alphafold_panels():
    """Fetches the model of each protein of the example pair, so they are cached beside the
    poster and the accession behind each screenshot is on record."""
    structures = strv.get_alphafold_structure(
        {protein: accession for _, _, protein, accession, _ in EXAMPLE_PAIR},
        output_dir=os.path.join(ASSET_DIR, 'structures'))

    for _, _, protein, accession, _ in EXAMPLE_PAIR:
        pdb_file, _, _, reason = structures[protein]
        state = os.path.basename(pdb_file) if pdb_file else f'unavailable ({reason})'
        print(f'  {protein} ({accession}): {state}')


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    config = utils.read_config(CONFIG_FILE)
    compare = load_app_page('app/pages/2_Compare_Parasites.py', 'poster_compare')

    print('1. interactions per parasite')
    build_interactions_per_parasite(config)

    print('1b. host proteins shared by many parasites')
    build_shared_proteins(config, compare)

    print('2. circos of shared host interactors')
    build_circos(config, compare)

    print(f'3. example interactome: {EXAMPLE_PARASITE} in {EXAMPLE_HOST}')
    example = example_predictions(config)
    build_body_figure(config, example)
    build_network(config, example)

    print('4. AlphaFold models')
    build_alphafold_panels()

    print(f'\nfigures in {OUT_DIR}')


if __name__ == '__main__':
    main()
