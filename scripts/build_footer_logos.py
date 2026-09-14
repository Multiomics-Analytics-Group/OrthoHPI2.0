'''
Prepare the source logos of the app footer (app/web_utils.py) at a common height, from the
copies in images/logos/source/.
'''
import os

import numpy as np
from PIL import Image, ImageChops

SOURCE_DIR = os.path.join('images', 'logos', 'source')
OUTPUT_DIR = os.path.join('images', 'logos')

# the height the footer draws them at; written at twice this for retina screens
LOGO_HEIGHT = 40

# colours kept when the logos are written as palette pngs
PALETTE_COLOURS = 256

# EggNOG's own dark teal, which its wordmark is reversed out of
EGGNOG_INK = (8, 45, 53)

# ebi.png holds the DeepMind and EBI logos side by side, with this gap between them
DEEPMIND_SPLIT = 414

# rows of tissues.png holding the wordmark; the tagline below dwarfs the other logos
TISSUES_WORDMARK = (16, 108)


def background_level(im):
    '''How light the background of a logo is, taken as the commonest colour in it.'''
    _, common = max(im.getcolors(im.width * im.height))

    return float(max(common))


def to_white_background(im):
    '''
    Lift a logo drawn on a light grey off it, by stretching the levels until the
    background is white.
    '''
    background = background_level(im)
    rgb = np.array(im, dtype=float) * 255.0 / background

    return Image.fromarray(np.clip(rgb, 0, 255).astype(np.uint8), 'RGB')


def invert_to_ink(im, ink):
    '''
    Turn a logo reversed out of a dark background into the same logo drawn in one colour
    on white, so that it can sit on the white footer beside the others.
    '''
    grey = np.array(im.convert('L'), dtype=float)

    # the background is the commonest colour rather than the darkest pixel
    _, common = max(im.getcolors(im.width * im.height))
    background = float(sum(common)) / 3.0
    weight = np.clip((grey - background) / max(255.0 - background, 1.0), 0, 1)[:, :, None]

    white = np.full(grey.shape + (3,), 255.0)
    rgb = white * (1.0 - weight) + np.array(ink, dtype=float) * weight

    return Image.fromarray(rgb.astype(np.uint8), 'RGB')


def trim(im):
    '''Drops the white margin around a logo, so the row is spaced by its own gap alone.'''
    white = Image.new('RGB', im.size, (255, 255, 255))
    box = ImageChops.difference(im, white).convert('L').getbbox()

    return im.crop(box) if box else im


def to_height(im, height):
    '''Scales a logo to a fixed height, keeping its aspect ratio.'''
    width = max(1, round(im.width * height / im.height))

    return im.resize((width, height), Image.LANCZOS)


def build():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    height = LOGO_HEIGHT * 2

    sources = {
        'eggnog': lambda im: invert_to_ink(im, EGGNOG_INK),
        # already on white
        'string': lambda im: im,
        'hpa': lambda im: im,
        'go': lambda im: im,
        'tissues': lambda im: to_white_background(
            im.crop((0, TISSUES_WORDMARK[0], im.width, TISSUES_WORDMARK[1]))),
        # the two logos of ebi.png are split so each is the size of the other logos
        'deepmind': lambda im: to_white_background(
            im.crop((0, 0, DEEPMIND_SPLIT, im.height))),
    }

    for name, prepare in sources.items():
        source = 'ebi' if name == 'deepmind' else name
        # composited onto white first, so antialiased edges are not blended against black
        original = Image.open(os.path.join(SOURCE_DIR, f'{source}.png')).convert('RGBA')
        flattened = Image.new('RGB', original.size, (255, 255, 255))
        flattened.paste(original, mask=original)

        out = to_height(trim(prepare(flattened)), height)
        # the footer inlines the files as data uris, and a palette halves the bytes
        out = out.quantize(colors=PALETTE_COLOURS)
        output_file = os.path.join(OUTPUT_DIR, f'{name}.png')
        out.save(output_file, optimize=True)
        print(f'{output_file}: {out.width}x{out.height}')


if __name__ == '__main__':
    build()
