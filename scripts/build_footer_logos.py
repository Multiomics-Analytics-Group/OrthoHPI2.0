"""
Prepare the source logos of the app footer (app/web_utils.py) at a common height, from the
copies in images/.

The originals were saved from each website at whatever size and background they happened to
be drawn on, and the footer drew them all at the same *width*. The aspect ratios run from
3.2:1 (EggNOG) to 10.7:1 (the Human Protein Atlas), so the Atlas rendered a third the
height of the rest, and the row showed EggNOG's dark navy plate and the grey plates of the
two screenshot crops beside two clean logos.

The logos are put on white rather than made transparent. Three of them (STRING, the Atlas,
Gene Ontology) are already drawn on white, and a grey `THE` of the Atlas wordmark cannot be
told apart from a black one at partial opacity, so keying them out guesses at the colour of
half the row. The footer is white, so matching it is the same thing to look at and nothing
has to be guessed. Only a dark theme would tell the difference, and the app does not offer
one.

Writes images/logos/<name>.png at LOGO_HEIGHT * 2 (retina), ready to be inlined by
app/web_utils.footer.

Run: .venv/bin/python scripts/build_footer_logos.py
"""
import os

import numpy as np
from PIL import Image, ImageChops

SOURCE_DIR = 'images'
OUTPUT_DIR = os.path.join('images', 'logos')

# the height the footer draws them at; the files are written at twice this so they stay
# sharp on a retina screen
LOGO_HEIGHT = 40

# colours kept when the logos are written as palette pngs
PALETTE_COLOURS = 256

# EggNOG's own dark teal, which its wordmark is reversed out of on the original
EGGNOG_INK = (8, 45, 53)

# the AlphaFold models are credited to both the people who built AlphaFold and the ones who
# host the database; ebi.png holds the two logos side by side, with this gap between them
DEEPMIND_SPLIT = 414

# rows of tissues.png holding the wordmark; the rest is the `Tissue expression database`
# tagline, which is set at its own scale and dwarfs the other logos once the whole crop is
# scaled to a common height
TISSUES_WORDMARK = (16, 108)


def background_level(im):
    """
    How light the background of a logo is, taken as the commonest colour in it. The
    screenshot crops are not on white -- the AlphaFold pair sits on a 243 grey and TISSUES
    on a 227 one -- so treating them as white would leave a visible plate behind them.

    :param im: RGB image
    :return: the lightest channel of the commonest colour, as a float
    """
    _, common = max(im.getcolors(im.width * im.height))

    return float(max(common))


def to_white_background(im):
    """
    Lift a logo drawn on a light grey off it, by stretching the levels until the background
    is white. The glyphs keep their colour and their antialiased edges fade into white
    rather than into the grey they were cut out with.

    :param im: RGB image
    :return: the same logo on white
    """
    background = background_level(im)
    rgb = np.array(im, dtype=float) * 255.0 / background

    return Image.fromarray(np.clip(rgb, 0, 255).astype(np.uint8), 'RGB')


def invert_to_ink(im, ink):
    """
    Turn a logo reversed out of a dark background into the same logo drawn in one colour on
    white, so that it can sit on the white footer beside the others.

    :param im: RGB image
    :param tuple ink: colour to draw the glyphs in
    :return: the recoloured logo
    """
    grey = np.array(im.convert('L'), dtype=float)

    # how far each pixel is from the background, where the glyphs are white. The background
    # is the commonest colour rather than the darkest pixel -- the wordmark has darker ink
    # inside its egg, and measuring from that leaves the plate as a pale wash of the ink.
    _, common = max(im.getcolors(im.width * im.height))
    background = float(sum(common)) / 3.0
    weight = np.clip((grey - background) / max(255.0 - background, 1.0), 0, 1)[:, :, None]

    white = np.full(grey.shape + (3,), 255.0)
    rgb = white * (1.0 - weight) + np.array(ink, dtype=float) * weight

    return Image.fromarray(rgb.astype(np.uint8), 'RGB')


def trim(im):
    """Drops the white margin around a logo, so the row is spaced by its own gap alone."""
    white = Image.new('RGB', im.size, (255, 255, 255))
    box = ImageChops.difference(im, white).convert('L').getbbox()

    return im.crop(box) if box else im


def to_height(im, height):
    """Scales a logo to a fixed height, keeping its aspect ratio."""
    width = max(1, round(im.width * height / im.height))

    return im.resize((width, height), Image.LANCZOS)


def build():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    height = LOGO_HEIGHT * 2

    sources = {
        'eggnog': lambda im: invert_to_ink(im, EGGNOG_INK),
        # already on white, and the three that need nothing done to them
        'string': lambda im: im,
        'hpa': lambda im: im,
        'go': lambda im: im,
        'tissues': lambda im: to_white_background(
            im.crop((0, TISSUES_WORDMARK[0], im.width, TISSUES_WORDMARK[1]))),
        # the two logos of ebi.png are drawn as one image, which left each of them at half
        # the size of every other logo in the row
        'deepmind': lambda im: to_white_background(
            im.crop((0, 0, DEEPMIND_SPLIT, im.height))),
    }

    for name, prepare in sources.items():
        source = 'ebi' if name == 'deepmind' else name
        # composited onto white first: string.png is transparent, and pasting it keeps the
        # antialiased edges from being blended against black
        original = Image.open(os.path.join(SOURCE_DIR, f'{source}.png')).convert('RGBA')
        flattened = Image.new('RGB', original.size, (255, 255, 255))
        flattened.paste(original, mask=original)

        out = to_height(trim(prepare(flattened)), height)
        # the footer inlines the files as data uris, so their size is paid on every page;
        # a wordmark has few enough colours that a palette is indistinguishable and half
        # the bytes
        out = out.quantize(colors=PALETTE_COLOURS)
        output_file = os.path.join(OUTPUT_DIR, f'{name}.png')
        out.save(output_file, optimize=True)
        print(f'{output_file}: {out.width}x{out.height}')


if __name__ == '__main__':
    build()
