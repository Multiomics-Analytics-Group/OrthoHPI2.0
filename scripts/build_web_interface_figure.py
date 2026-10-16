'''
The web interface figure: screenshots of the network page laid out in two columns, each
column's screenshots at one width and the columns as tall as each other, lettered down
the first column and then the second.

  A  the body figure and the query controls
  B  the network
  C  the AlphaFold models of a link's proteins, opened by clicking it
  D  the GO biological processes enriched among the host proteins

The screenshots are paper/figures/web/*.png, taken from the running app by
scripts/capture_web_screenshots.py. Writes paper/figures/web_interface.svg/.png.
'''
import argparse
import os
import string
import sys

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts import figure_style

# the screenshots, column by column
COLUMNS = [['controls', 'network'], ['alphafold', 'enrichment']]
# the name each screenshot's panel is titled with, after its letter
TITLES = {
    'controls': 'Query',
    'network': 'Network',
    'alphafold': 'AlphaFold models of a link',
    'enrichment': 'Enriched GO processes',
}
# inches between the columns, between panels down a column, and for the title above each
GAP = 0.2
ROW_GAP = 0.08
TITLE_ROOM = 0.2
BORDER_COLOR = '#c7cfd9'


def column_widths(images):
    '''The width of each column that makes the columns equally tall across the figure width.'''
    # the screenshot height of each column per inch of its width
    tallness = [sum(image.shape[0] / image.shape[1] for image in column) for column in images]
    room = figure_style.WIDTH - GAP * (len(images) - 1)
    # widths inversely proportional to tallness, summing to the room
    inverse = [1 / t for t in tallness]

    return [room * i / sum(inverse) for i in inverse]


def draw(columns, screenshot_dir, output_stem):
    figure_style.apply()
    images = [[mpimg.imread(os.path.join(screenshot_dir, f'{name}.png')) for name in column]
              for column in columns]
    widths = column_widths(images)
    panels = max(len(column) for column in columns)
    chrome = (TITLE_ROOM + ROW_GAP) * panels - ROW_GAP
    height = widths[0] * sum(image.shape[0] / image.shape[1] for image in images[0])
    # a figure taller than a page shrinks every column alike, leaving them narrower
    scale = min(1, (figure_style.MAX_HEIGHT - chrome) / height)
    widths = [w * scale for w in widths]
    fig_height = height * scale + chrome
    fig_width = sum(widths) + GAP * (len(widths) - 1)
    fig = plt.figure(figsize=(fig_width, fig_height))

    letters = iter(string.ascii_uppercase)
    x = 0
    for column, column_images, width in zip(columns, images, widths):
        top = fig_height
        for name, image in zip(column, column_images):
            top -= TITLE_ROOM
            panel_height = width * image.shape[0] / image.shape[1]
            ax = fig.add_axes([x / fig_width, (top - panel_height) / fig_height,
                               width / fig_width, panel_height / fig_height])
            ax.imshow(image, interpolation='lanczos')
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_color(BORDER_COLOR)
                spine.set_linewidth(0.6)
            ax.set_title(f'{next(letters)}  {TITLES[name]}', loc='left', fontweight='bold', pad=3)
            top -= panel_height + ROW_GAP
        x += width + GAP

    for extension in ('svg', 'png'):
        fig.savefig(f'{output_stem}.{extension}', dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--screenshots', default='paper/figures/web')
    parser.add_argument('--figures', default='paper/figures', help='output directory')
    args = parser.parse_args()

    output_stem = os.path.join(args.figures, 'web_interface')
    draw(COLUMNS, args.screenshots, output_stem)
    print(f"Wrote {output_stem}.svg and .png")
