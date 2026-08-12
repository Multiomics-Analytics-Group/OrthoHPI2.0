"""
Assembles poster/poster-v2.svg, the A0 poster, from the figures build_poster_figures.py
writes and the logos in poster/assets/.

Written as a generator rather than by hand so the poster can be rebuilt when a figure
changes, and so the headline numbers are read out of the data instead of being typed in
and quietly going stale.

The document is A0 portrait with a viewBox of 0 0 841 1189, so one user unit is one
millimetre and every size in this file is a size on the printed sheet. Figures are inlined
as groups rather than referenced with <image>, which keeps them vector and editable and
avoids Inkscape's patchy handling of nested SVG; their internal ids are prefixed on the
way in so two figures cannot collide over a clipPath. Rasters are embedded as data URIs so
the file stands alone.

The output is meant to be opened in Inkscape: the bands are Inkscape layers, and all text
is real text.

Run build_poster_figures.py first, then this. The print file is exported from the result:

    ~/.venvs/streamlit/bin/python scripts/build_poster_figures.py
    ~/.venvs/streamlit/bin/python scripts/build_poster.py
    /Applications/Inkscape.app/Contents/MacOS/inkscape --export-type=pdf \\
        --export-text-to-path --export-filename=poster/poster-v2.pdf poster/poster-v2.svg

--export-text-to-path is what stops a print shop substituting a font it does not have.
Rebuilding overwrites poster/poster-v2.svg, so anything edited into it by hand in Inkscape
is lost -- port such edits back into this file instead.
"""
import base64
import os
import re
import sys
import textwrap
import xml.etree.ElementTree as ET

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'app'))

import logging  # noqa: E402
logging.getLogger('streamlit').setLevel(logging.ERROR)

import utils  # noqa: E402
import web_utils  # noqa: E402

FIG_DIR = os.path.join(REPO, 'poster', 'figures')
ASSET_DIR = os.path.join(REPO, 'poster', 'assets')
IMAGE_DIR = os.path.join(REPO, 'images')
OUT_FILE = os.path.join(REPO, 'poster', 'poster-v2.svg')

PAGE_W, PAGE_H = 841, 1189
MARGIN = 25
CONTENT_W = PAGE_W - 2 * MARGIN          # 791
RIGHT = MARGIN + CONTENT_W               # 816

# Neutral chrome. The figures carry the Okabe-Ito clade colours and the blues of the chord
# ramp; if the poster's own furniture were coloured too, those would stop reading as data.
# The one accent is the dark blue the app already sets its headings in.
INK = '#1d2126'
BODY = '#3c434c'
MUTED = '#6b7280'
RULE = '#ccd2d8'
ACCENT = '#023858'
HEADER_BG = '#22272e'      # the DTU and Novo Nordisk marks are white, so the band is dark
PANEL = '#f4f6f8'

FONT = 'Arial, Helvetica, sans-serif'

# Vertical plan. Each band is (top, height); the gaps between them are what is left over.
BANDS = {
    'header': (0, 132),
    'intro': (146, 108),
    'pipeline': (266, 172),
    'figures': (450, 658),
    'footer': (1120, 60),
}


# ----------------------------------------------------------------------------------
# primitives
# ----------------------------------------------------------------------------------

def escape(text):
    return (str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;'))


def text(x, y, content, size=5, fill=BODY, weight='normal', anchor='start',
         style='normal', spacing=None, family=FONT):
    attrs = [f'x="{x:.2f}"', f'y="{y:.2f}"', f'font-family="{family}"',
             f'font-size="{size}"', f'fill="{fill}"']
    if weight != 'normal':
        attrs.append(f'font-weight="{weight}"')
    if anchor != 'start':
        attrs.append(f'text-anchor="{anchor}"')
    if style != 'normal':
        attrs.append(f'font-style="{style}"')
    if spacing:
        attrs.append(f'letter-spacing="{spacing}"')

    return f'<text {" ".join(attrs)}>{escape(content)}</text>\n'


def paragraph(x, y, content, width, size=5, fill=BODY, leading=1.35, weight='normal'):
    """
    Text wrapped to a column, as a run of <tspan> lines.

    The wrap point is worked out from an average glyph advance rather than measured -- there
    is no font engine here -- which is close enough for Arial at these sizes and leaves the
    lines as editable text instead of baking them into paths.
    """
    per_character = size * 0.52
    columns = max(int(width / per_character), 8)
    lines = textwrap.wrap(content, width=columns) or ['']

    spans = ''.join(
        f'<tspan x="{x:.2f}" dy="{0 if i == 0 else size * leading:.2f}">{escape(line)}</tspan>'
        for i, line in enumerate(lines))
    weight_attr = f' font-weight="{weight}"' if weight != 'normal' else ''

    return (f'<text x="{x:.2f}" y="{y:.2f}" font-family="{FONT}" font-size="{size}" '
            f'fill="{fill}"{weight_attr}>{spans}</text>\n'), len(lines) * size * leading


def bullet(x, y, content, width, size=5, fill=BODY):
    """A bullet and its text, the text indented so its wrapped lines line up under
    themselves rather than under the bullet."""
    indent = size * 1.1
    dot = (f'<circle cx="{x + size * 0.35:.2f}" cy="{y - size * 0.32:.2f}" '
           f'r="{size * 0.16:.2f}" fill="{MUTED}"/>\n')
    body, height = paragraph(x + indent, y, content, width - indent, size=size, fill=fill)

    return dot + body, height


def rounded_panel(x, y, width, height, fill=PANEL, stroke='none', radius=3, dash=None):
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ''
    stroke_attr = f' stroke="{stroke}" stroke-width="0.6"{dash_attr}' if stroke != 'none' else ''

    return (f'<rect x="{x:.2f}" y="{y:.2f}" width="{width:.2f}" height="{height:.2f}" '
            f'rx="{radius}" fill="{fill}"{stroke_attr}/>\n')


def section_heading(y, label):
    """A section title over a hairline. v1 used a filled bar the width of the poster for
    each of these; three of them cost most of a panel's worth of height and made the page
    look like a stack of headers rather than a stack of content."""
    out = text(MARGIN, y, label.upper(), size=15, fill=ACCENT, weight='bold', spacing='2')
    out += (f'<line x1="{MARGIN}" y1="{y + 6:.2f}" x2="{RIGHT}" y2="{y + 6:.2f}" '
            f'stroke="{ACCENT}" stroke-width="0.9"/>\n')

    return out


def caption(x, y, label, body_text, width, number=None):
    """The heading over a figure: an optional step number, the thing the figure shows, and
    a sentence on how to read it."""
    out = ''
    offset = 0
    if number is not None:
        out += (f'<circle cx="{x + 5:.2f}" cy="{y - 2.6:.2f}" r="5" fill="{ACCENT}"/>\n')
        out += text(x + 5, y, str(number), size=7, fill='white', weight='bold',
                    anchor='middle')
        offset = 13.5

    out += text(x + offset, y, label, size=9.5, fill=INK, weight='bold')
    body, _ = paragraph(x, y + 11.5, body_text, width, size=6.4, fill=MUTED)

    return out + body


# ----------------------------------------------------------------------------------
# bringing figures and images in
# ----------------------------------------------------------------------------------

SVG_NS = 'http://www.w3.org/2000/svg'
XLINK_NS = 'http://www.w3.org/1999/xlink'
ET.register_namespace('', SVG_NS)
ET.register_namespace('xlink', XLINK_NS)


def namespace_ids(markup, prefix):
    """
    Prefixes every id in a figure and every reference to one.

    Plotly names its clip paths by figure index, so two of its exports side by side would
    both define `clip1` and the second would silently win, cropping the first to the wrong
    rectangle. Renaming on the way in is what makes several figures able to share a
    document at all.
    """
    ids = set(re.findall(r'\bid="([^"]+)"', markup))
    for original in sorted(ids, key=len, reverse=True):
        renamed = f'{prefix}-{original}'
        markup = markup.replace(f'id="{original}"', f'id="{renamed}"')
        markup = markup.replace(f'url(#{original})', f'url(#{renamed})')
        markup = markup.replace(f'href="#{original}"', f'href="#{renamed}"')

    return markup


def inline_figure(path, x, y, width, height, prefix, align='middle'):
    """
    Places a figure inside a box, scaled to fit and keeping its aspect.

    The figure's own viewBox is what it is measured by, so a figure written in millimetres
    lands at exactly the size it was written for and one carrying foreign coordinates --
    Plotly's pixels, the body drawing's user units -- is scaled into the same box.
    """
    tree = ET.parse(path)
    root = tree.getroot()
    view_box = root.get('viewBox')
    if view_box:
        vx, vy, vw, vh = (float(v) for v in view_box.replace(',', ' ').split())
    else:
        vx, vy = 0.0, 0.0
        vw, vh = float(root.get('width', width)), float(root.get('height', height))

    scale = min(width / vw, height / vh)
    dx = x + (width - vw * scale) / 2 - vx * scale
    dy = y + (0 if align == 'top' else (height - vh * scale) / 2) - vy * scale

    inner = ''.join(ET.tostring(child, encoding='unicode') for child in root)
    inner = namespace_ids(inner, prefix)
    # ElementTree writes the default namespace onto every child it serialises
    inner = inner.replace(f'xmlns:ns0="{SVG_NS}"', '').replace('ns0:', '')

    return (f'<g id="{prefix}" transform="translate({dx:.3f},{dy:.3f}) '
            f'scale({scale:.5f})">\n{inner}\n</g>\n')


def embed_image(path, x, y, width, height, preserve='xMidYMid meet'):
    """A raster as a data URI, so the poster is one self-contained file."""
    with open(path, 'rb') as handle:
        payload = base64.b64encode(handle.read()).decode('ascii')
    kind = 'jpeg' if path.lower().endswith(('.jpg', '.jpeg')) else 'png'

    return (f'<image x="{x:.2f}" y="{y:.2f}" width="{width:.2f}" height="{height:.2f}" '
            f'preserveAspectRatio="{preserve}" '
            f'xlink:href="data:image/{kind};base64,{payload}"/>\n')


def image_aspect(path):
    from PIL import Image
    with Image.open(path) as image:
        return image.width / image.height


# ----------------------------------------------------------------------------------
# the bands
# ----------------------------------------------------------------------------------

# spelled as v1 spells it, which is the version that was called final
TITLE = 'OrthoHPI2.0: Prediction of host-parasite protein-protein interactions'
MONA = 'Multi-omics Network Analytics group'
AUTHORS = ('Ida K.S. Meitil¹, Yesid Cuesta Astroz², '
           'Juan Esteban Vargas Gomez², Alberto Santos¹')
AFFILIATIONS = [
    '1  Department of health technology, Technical University of Denmark, Kgs. Lyngby, Denmark',
    '2  Escuela de microbiología, Universidad de Antioquia, Medellin, Columbia',
]


def build_header():
    top, height = BANDS['header']
    out = f'<rect x="0" y="{top}" width="{PAGE_W}" height="{height}" fill="{HEADER_BG}"/>\n'

    # the two funder marks on the left, the MONA mark on the right, the title between them
    dtu = os.path.join(IMAGE_DIR, 'dtu-logo-white.png')
    nnf = os.path.join(IMAGE_DIR, 'nnf-logo-white.png')
    mona = os.path.join(ASSET_DIR, 'mona-logo.png')
    out += embed_image(dtu, MARGIN, top + 16, 56, 39, preserve='xMinYMid meet')
    out += embed_image(nnf, MARGIN, top + 60, 56, 51, preserve='xMinYMid meet')
    out += embed_image(mona, PAGE_W - MARGIN - 116, top + 10, 116, 80, preserve='xMaxYMid meet')
    # the group's name, which in v1 was live text beside the mark rather than part of it and
    # so did not survive being lifted out of the PDF
    out += text(PAGE_W - MARGIN, top + 106, MONA, size=6.8, fill='#c8ced6', anchor='end',
                spacing='0.4')

    left = MARGIN + 70
    # the title block hangs lower than the logos so the header reads as optically centred
    title, used = paragraph(left, top + 45.35, TITLE, 580, size=30, fill='white',
                            weight='bold', leading=1.14)
    out += title
    out += text(left, top + 106.4, AUTHORS, size=9.5, fill='#c8ced6')
    for i, line in enumerate(AFFILIATIONS):
        out += text(left, top + 119.4 + i * 9, line, size=7.4, fill='#98a2ad',
                    style='italic')

    return out


INTRO_MAIN = [
    'OrthoHPI predicts host-parasite protein-protein interactions (PPIs) via '
    'orthology-transfer of intra-species PPIs to inter-species PPIs.',
    'Filters based on host protein expression data and protein localization',
]
INTRO_NEW = [
    'Four hosts: human, pig, mouse and rat',
    '40 parasites',
    'Cell-type level expression data',
    'DeepLoc2.1 [1] for protein localization prediction',
]


def build_intro(stats):
    top, height = BANDS['intro']
    out = section_heading(top + 9, 'Introduction')

    y = top + 32
    column = 250
    for item in INTRO_MAIN:
        block, used = bullet(MARGIN, y, item, column, size=8.5)
        out += block
        y += used + 2.5

    x2 = MARGIN + column + 20
    out += text(x2, top + 32, 'New in this version:', size=8.5, fill=INK, weight='bold')
    y = top + 46
    for item in INTRO_NEW:
        block, used = bullet(x2, y, item, 205, size=8.5)
        out += block
        y += used + 1.5

    # the headline numbers, which v1 left buried in the bullets above
    x3 = MARGIN + column + 20 + 205 + 24
    tiles = [(f"{stats['parasites']}", 'parasites'),
             (f"{stats['interactions']:,}", 'predicted interactions')]
    tile_w = (RIGHT - x3 - 2 * 6) / 3
    for i, (value, label) in enumerate(tiles):
        tx = x3 + i * (tile_w + 6)
        out += rounded_panel(tx, top + 24, tile_w, 66)
        out += text(tx + tile_w / 2, top + 58, value, size=26, fill=ACCENT, weight='bold',
                    anchor='middle')
        body, _ = paragraph(tx + 4, top + 72, label, tile_w - 8, size=6.4, fill=MUTED)
        # centred under the number
        out += body.replace(f'x="{tx + 4:.2f}"', f'x="{tx + tile_w / 2:.2f}"').replace(
            'font-size="6.4"', 'font-size="6.4" text-anchor="middle"')

    # the third tile is the host split, which the bar chart used to carry as three columns
    tx = x3 + 2 * (tile_w + 6)
    out += rounded_panel(tx, top + 24, tile_w, 66)
    out += text(tx + 7, top + 36, 'by host', size=6, fill=MUTED)
    out += text(tx + tile_w - 7, top + 36, 'parasites · PPIs', size=5.2, fill=MUTED,
                anchor='end')
    for i, (host, parasites, interactions) in enumerate(stats['by_host']):
        row_y = top + 48 + i * 12
        out += text(tx + 7, row_y, host, size=6.4, fill=INK, style='italic')
        out += text(tx + tile_w - 7, row_y, f'{parasites} · {interactions:,}', size=6.4,
                    fill=ACCENT, weight='bold', anchor='end')

    return out


def build_pipeline():
    """
    The four stages of the prediction, left to right.

    v1 set these as a 2x2 with the filter hanging off to one side, which left a band of
    empty sheet under it and made the reader work out that the filter came third. As one
    row the order is the order they are read in, and the band it used to waste goes to the
    figures.
    """
    top, height = BANDS['pipeline']
    out = section_heading(top + 9, 'Prediction pipeline')

    box_y = top + 30
    box_h = height - 36
    stage_w = 176
    gap = (CONTENT_W - 4 * stage_w) / 3      # 29mm, the arrows live here

    stages = [
        ('Step 1', 'Known interaction', 'An interaction measured within one species'),
        ('Step 2', 'Orthology transfer', 'Both proteins sit in orthologous groups'),
        ('Filter', 'Which proteins qualify', None),
        ('Step 3', 'Host-parasite interaction', 'The parasite and host members inherit the link'),
    ]

    for i, (step, title, subtitle) in enumerate(stages):
        x = MARGIN + i * (stage_w + gap)
        dashed = step == 'Filter'
        out += rounded_panel(x, box_y, stage_w, box_h, fill='white' if dashed else PANEL,
                             stroke='#1b7837' if dashed else 'none',
                             dash='3 2.4' if dashed else None)
        out += text(x + 10, box_y + 16, f'{step.upper()}  ·  {title.upper()}', size=6.4,
                    fill='#1b7837' if dashed else MUTED, weight='bold', spacing='0.5')
        if subtitle:
            # 6.4 rather than 6.8 so the longest of the four -- step 3's -- still sets on
            # one line; one stage wrapping while the others do not makes the row look ragged
            body, _ = paragraph(x + 10, box_y + 30, subtitle, stage_w - 20, size=8.5, fill=INK,
                                weight='bold')
            out += body

        if i < len(stages) - 1:
            ax = x + stage_w + 4
            ay = box_y + box_h / 2
            out += (f'<path d="M{ax:.2f},{ay:.2f} L{ax + gap - 8:.2f},{ay:.2f}" '
                    f'stroke="{MUTED}" stroke-width="1.1" marker-end="url(#arrow)"/>\n')

        out += stage_art(i, x, box_y, stage_w, box_h)

    return out


def stage_art(index, x, y, width, height):
    """The little diagram inside each stage, which is what the stage actually says. The
    wording is v1's; only the arrangement changed."""
    cx = x + width / 2
    out = ''

    if index == 0:
        art_y = y + 80
        out += text(cx, art_y - 21, 'interact in STRING 12.0', size=6.4, fill=MUTED,
                    anchor='middle')
        for dx, label in [(-30, 'A'), (30, 'B')]:
            out += (f'<circle cx="{cx + dx:.2f}" cy="{art_y:.2f}" r="11" fill="white" '
                    f'stroke="{INK}" stroke-width="1.2"/>\n')
            out += text(cx + dx, art_y + 4.4, label, size=12, fill=INK, weight='bold',
                        anchor='middle')
        out += (f'<line x1="{cx - 19:.2f}" y1="{art_y:.2f}" x2="{cx + 19:.2f}" '
                f'y2="{art_y:.2f}" stroke="{INK}" stroke-width="1.6"/>\n')
        out += text(cx, art_y + 27, 'experimental or database evidence', size=6.4, fill=MUTED,
                    anchor='middle')

    elif index == 1:
        art_y = y + 90
        body, _ = paragraph(x + 10, y + 56, 'EggNOG 6 assigns every protein to an '
                            'orthologous group', width - 20, size=6.4, fill=MUTED)
        out += body
        for dx, label, sub in [(-44, 'Group 1', 'orthologs of A'), (44, 'Group 2', 'orthologs of B')]:
            out += rounded_panel(cx + dx - 38, art_y - 15, 76, 30, fill='#eaf2fa',
                                 stroke='#1f6fb2')
            out += text(cx + dx, art_y - 2, label, size=7.4, fill='#1f6fb2', weight='bold',
                        anchor='middle')
            out += text(cx + dx, art_y + 9, sub, size=6.4, fill=MUTED, anchor='middle')
        out += (f'<line x1="{cx - 8:.2f}" y1="{art_y:.2f}" x2="{cx + 8:.2f}" y2="{art_y:.2f}" '
                f'stroke="#1f6fb2" stroke-width="1.6"/>\n')
        out += text(cx, art_y + 29, 'group–group link', size=6.4, fill=MUTED, anchor='middle')
        out += text(cx, art_y + 41, 'score ≥ 0.7', size=6.6, fill=INK, anchor='middle',
                    family='Courier New, monospace')

    elif index == 2:
        # the two criteria sit where the diagrams of the other three stages do, so the row
        # of stages reads as one band rather than one box that is top-heavy
        art_y = y + 50
        for dx, heading, colour, detail in [
                (-42, 'Parasite protein', '#D55E00',
                 'predicted secreted or membrane-bound'),
                (42, 'Host protein', INK,
                 'Predicted surface-exposed, expressed in a tissue the parasite infects')]:
            out += text(cx + dx, art_y + 12, heading, size=7.4, fill=colour, weight='bold',
                        anchor='middle')
            body, _ = paragraph(cx + dx - 42, art_y + 24, detail, 84, size=6.4, fill=BODY)
            out += body.replace(f'x="{cx + dx - 42:.2f}"', f'x="{cx + dx:.2f}"').replace(
                'font-size="6.4"', 'font-size="6.4" text-anchor="middle"')

    elif index == 3:
        art_y = y + 80
        out += text(cx, art_y - 21, 'predicted interaction', size=6.4, fill='#1b7837',
                    anchor='middle')
        out += (f'<rect x="{cx - 42:.2f}" y="{art_y - 10:.2f}" width="20" height="20" '
                f'fill="#D55E00" transform="rotate(45 {cx - 32:.2f} {art_y:.2f})"/>\n')
        out += text(cx - 32, art_y + 4.4, 'P', size=11, fill='white', weight='bold',
                    anchor='middle')
        out += (f'<circle cx="{cx + 32:.2f}" cy="{art_y:.2f}" r="11" fill="#3f4650"/>\n')
        out += text(cx + 32, art_y + 4.4, 'H', size=11, fill='white', weight='bold',
                    anchor='middle')
        out += (f'<line x1="{cx - 18:.2f}" y1="{art_y:.2f}" x2="{cx + 19:.2f}" y2="{art_y:.2f}" '
                f'stroke="#1b7837" stroke-width="1.6" stroke-dasharray="4 2.6"/>\n')
        for dx, label, sub in [(-32, 'parasite protein', 'member of Group 1'),
                               (32, 'host protein', 'member of Group 2')]:
            out += text(cx + dx, art_y + 27, label, size=6.4, fill=BODY, anchor='middle')
            out += text(cx + dx, art_y + 36, sub, size=6, fill=MUTED, anchor='middle')

    return out


def build_figures():
    """The web resource, as a zoom: the whole resource, then the parasites compared, then
    one interactome, then one interacting pair."""
    top, band_height = BANDS['figures']
    out = section_heading(top + 9, 'The web resource')

    # The three rows and what each spends above its figures on a caption. Kept as data so
    # the band can be checked against the space it has: the first version of this ran the
    # last row 19mm into the footer, which is not the kind of thing to find on a printed A0.
    ROWS = [(22, 155), (22, 290), (30, 95)]
    GAPS = [12, 10]
    needed = 22 + sum(c + h for c, h in ROWS) + sum(GAPS)
    if needed > band_height:
        raise SystemExit(f'figures need {needed}mm but the band is {band_height}mm; '
                         f'adjust ROWS, GAPS or BANDS["figures"]')

    # 1 -- the scale of the resource
    y = top + 22
    caption_h, figure_h = ROWS[0]
    out += caption(MARGIN, y, 'Host proteins many parasites reach',
                   'The 12 most shared human proteins. Dot size is how many of that '
                   'parasite\u2019s proteins reach it.', 600, number=1)
    out += inline_figure(os.path.join(FIG_DIR, 'fig1b-shared-proteins.svg'),
                         MARGIN, y + caption_h, 780, figure_h, 'fig1b')

    # 2 and 3 -- the parasites compared, and one interactome opened up
    y = y + caption_h + figure_h + GAPS[0]
    caption_h, figure_h = ROWS[1]
    out += caption(MARGIN, y, 'Which parasites reach the same host proteins',
                   'Pairs sharing at least 10 host targets, in human. Darker chords share more.',
                   370, number=2)
    out += inline_figure(os.path.join(FIG_DIR, 'fig2-circos.svg'),
                         MARGIN, y + caption_h, 380, figure_h, 'fig2')

    x_net = MARGIN + 380 + 16
    out += caption(x_net, y, 'One interactome',
                   'P. falciparum against human at confidence ≥ 0.6, not filtered by tissue.',
                   395, number=3)
    out += inline_figure(os.path.join(FIG_DIR, 'fig3b-network.svg'),
                         x_net, y + caption_h, 395, figure_h, 'fig3b')

    # 4 and 5 -- where that interactome acts, and one of its pairs down to structure
    y = y + caption_h + figure_h + GAPS[1]
    caption_h, figure_h = ROWS[2]
    out += caption(MARGIN, y, 'Where it acts',
                   'The organs P. falciparum infects, shaded by interactions reaching them.',
                   170, number=4)
    out += inline_figure(os.path.join(FIG_DIR, 'fig3a-body.svg'),
                         MARGIN, y + caption_h, 170, figure_h, 'fig3a')

    x_af = MARGIN + 170 + 16
    out += caption(x_af, y, 'Structures of an interacting pair',
                   'AlphaFold models of one predicted interaction.',
                   580, number=5)
    out += alphafold_panel(x_af, y + caption_h, RIGHT - x_af, figure_h)

    return out


def alphafold_panel(x, y, width, height):
    """
    The two AlphaFold models side by side, each labelled with the protein it is.

    The models are screenshots of the app's 3Dmol viewer, saved into images/. The caption
    sits beside the model rather than over it: stacked, the text took a third of a panel
    only 95mm tall and left the fold too small to make out.
    """
    out = ''
    panels = [
        ('PF3D7-AF.png', 'PARASITE', '#D55E00', 'Plasmodium falciparum',
         'PF3D7_1104100', 'Synaptobrevin, putative · Q8IIW1'),
        ('YKT6-AF.png', 'HOST', '#3f4650', 'Homo sapiens',
         'YKT6', 'Synaptobrevin homolog YKT6 · O15498'),
    ]
    half = (width - 14) / 2
    caption_w = 104
    for i, (filename, tag, colour, species, protein, description) in enumerate(panels):
        px = x + i * (half + 14)
        out += rounded_panel(px, y, half, height, fill=PANEL)

        # the pill follows its word: at a fixed width HOST sat in a puddle and PARASITE
        # ran out of the end of it
        badge_size = 5.2
        badge_w = len(tag) * badge_size * 0.68 + 7
        out += (f'<rect x="{px + 7:.2f}" y="{y + 8:.2f}" width="{badge_w:.2f}" height="7" '
                f'rx="3.5" fill="{colour}"/>\n')
        out += text(px + 7 + badge_w / 2, y + 13, tag, size=badge_size, fill='white',
                    weight='bold', anchor='middle', spacing='0.3')
        out += text(px + 8, y + 30, species, size=7, fill=INK, weight='bold', style='italic')
        out += text(px + 8, y + 41, protein, size=8, fill=INK, weight='bold')
        body, _ = paragraph(px + 8, y + 52, description, caption_w - 12, size=6.4, fill=MUTED)
        out += body

        out += embed_image(os.path.join(IMAGE_DIR, filename), px + caption_w, y + 4,
                           half - caption_w - 6, height - 8)

    return out


def build_footer():
    top, height = BANDS['footer']
    out = (f'<line x1="{MARGIN}" y1="{top}" x2="{RIGHT}" y2="{top}" '
           f'stroke="{RULE}" stroke-width="0.8"/>\n')

    out += text(MARGIN, top + 13, 'Developed with data from:', size=6.4, fill=MUTED)

    logos = ['string.png', 'eggnog.png', 'go.png', 'tissues.png', 'hpa.png', 'ebi.png']
    lx = MARGIN
    ly = top + 18
    for filename in logos:
        path = os.path.join(IMAGE_DIR, filename)
        if not os.path.exists(path):
            continue
        logo_h = 13
        logo_w = min(logo_h * image_aspect(path), 74)
        out += embed_image(path, lx, ly, logo_w, logo_h, preserve='xMinYMid meet')
        lx += logo_w + 9

    # the two codes, right aligned. Placed before the references so their left edge is
    # what the reference column is allowed to run up to
    qr = [('qr-github.png', 'Code'), ('qr-streamlit.png', 'App')]
    qx = RIGHT - len(qr) * 52

    out += text(MARGIN, top + 50, 'References:', size=6.4, fill=MUTED)
    ry = top + 50
    for i, reference in enumerate(REFERENCES, start=1):
        entry, used = paragraph(MARGIN + 40, ry, f'[{i}] {reference}', qx - MARGIN - 48,
                                size=6, fill=MUTED)
        out += entry
        ry += used + 1.5
    if ry > top + height:
        raise SystemExit(f'references run {ry - (top + height):.0f}mm past the footer; '
                         f'grow BANDS["footer"] and move it up')

    for filename, label in qr:
        path = os.path.join(IMAGE_DIR, filename)
        if not os.path.exists(path):
            continue
        out += embed_image(path, qx, top + 5, 46, 46)
        out += text(qx + 23, top + 57, label, size=6.4, fill=MUTED, anchor='middle')
        qx += 52

    return out


# Numbered in the order they are listed. The footer grows downwards from a fixed top, so
# build_footer stops the build rather than letting a long list run off the sheet.
REFERENCES = [
    'Ødum MT, Teufel F, Thumuluri V, et al. DeepLoc 2.1: multi-label membrane protein '
    'type prediction using protein language models. Nucleic Acids Res. '
    '2024;52(W1):W215–W220. doi:10.1093/nar/gkae237',
]


# ----------------------------------------------------------------------------------

def headline_stats():
    """The three numbers on the poster, read out of the configuration and the predictions
    so they cannot drift from what was actually built."""
    config = utils.read_config(os.path.join(REPO, 'config.yml'))
    predictions = web_utils.load_predictions(os.path.join(REPO, 'data'))

    by_host = []
    for host, taxids in web_utils.get_host_groups(config, predictions).items():
        rows = predictions[predictions['taxid2'].isin(taxids)]
        by_host.append((host, rows['taxid1_label'].nunique(),
                        len(rows.drop_duplicates(subset=['source', 'target']))))

    return {
        'hosts': len(config['hosts']),
        'parasites': predictions['taxid1_label'].nunique(),
        'interactions': len(predictions.drop_duplicates(subset=['source', 'target'])),
        'by_host': by_host,
    }


LAYERS = [
    ('header', build_header),
    ('intro', build_intro),
    ('pipeline', build_pipeline),
    ('figures', build_figures),
    ('footer', build_footer),
]


def main():
    stats = headline_stats()
    print(f'  {stats["hosts"]} hosts, {stats["parasites"]} parasites, '
          f'{stats["interactions"]:,} predicted interactions')

    parts = [
        f'<svg xmlns="{SVG_NS}" xmlns:xlink="{XLINK_NS}" '
        f'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
        f'width="{PAGE_W}mm" height="{PAGE_H}mm" viewBox="0 0 {PAGE_W} {PAGE_H}">\n',
        '<defs>\n'
        '<marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5.5" '
        f'markerHeight="5.5" orient="auto-start-reverse">\n'
        f'<path d="M0,1 L10,5 L0,9 z" fill="{MUTED}"/>\n</marker>\n</defs>\n',
        f'<rect width="{PAGE_W}" height="{PAGE_H}" fill="white"/>\n',
    ]

    for name, builder in LAYERS:
        content = builder(stats) if name == 'intro' else builder()
        parts.append(f'<g inkscape:groupmode="layer" inkscape:label="{name}" id="layer-{name}">\n')
        parts.append(content)
        parts.append('</g>\n')

    parts.append('</svg>\n')

    os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
    with open(OUT_FILE, 'w') as handle:
        handle.write(''.join(parts))

    print(f'  wrote {os.path.relpath(OUT_FILE, REPO)} '
          f'({os.path.getsize(OUT_FILE) / 1024 / 1024:.1f} MB)')


if __name__ == '__main__':
    main()
