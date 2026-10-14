'''
The two web interface figures: screenshots of the app laid out in rows, each row as tall
as fits the page width with its screenshots at one height, lettered in reading order.

  web_overview  the app across parasites and hosts
    A  the home page: the hosts, and the parasites by taxonomic group
    B  the multi-host parasites page: interactions shared between the hosts of a parasite
    C  the parasites-of-a-host page: host interactors shared by each pair of parasites,
       and the dialog listing them for one pair
  web_network   one host-parasite pair on the network page
    A  the query controls and the network
    B  the AlphaFold models of a link's proteins, opened by clicking it
    C  the cell types expressing the host proteins
    D  the GO biological processes enriched in the network

The screenshots are paper/figures/web/*.png, taken from the running app by
scripts/capture_web_screenshots.py. Writes paper/figures/web_overview.svg/.png and
paper/figures/web_network.svg/.png.
'''
import argparse
import os
import string
import sys

import matplotlib.image as mpimg
from matplotlib.patches import ConnectionPatch
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts import figure_style

# the screenshots of each figure, row by row
FIGURES = {
    'web_overview': [['home', 'multi_host'], ['parasite_overlap']],
    'web_network': [['network', 'alphafold'], ['cell_types', 'enrichment']],
}
# a click in the first screenshot opens the second
ARROWS = {'web_network': [('network', 'alphafold')]}
# the name each screenshot's panel is titled with, after its letter
TITLES = {
    'home': 'Home page',
    'multi_host': 'Multi-host parasites page',
    'parasite_overlap': 'Parasites of a host page',
    'network': 'Network page',
    'alphafold': 'AlphaFold models of a link',
    'cell_types': 'Host cell types',
    'enrichment': 'Enriched GO processes',
}
# inches between panels across a row, and for the letter above each screenshot
GAP = 0.2
TITLE_ROOM = 0.2
BORDER_COLOR = '#c7cfd9'


def row_heights(images):
    '''The screenshot height of each row that makes it fill the figure width.'''
    heights = []
    for row in images:
        aspects = [image.shape[1] / image.shape[0] for image in row]
        heights.append((figure_style.WIDTH - GAP * (len(row) - 1)) / sum(aspects))

    return heights


def draw(rows, arrows, screenshot_dir, output_stem):
    figure_style.apply()
    images = [[mpimg.imread(os.path.join(screenshot_dir, f'{name}.png')) for name in row]
              for row in rows]
    heights = row_heights(images)
    # a figure taller than a page shrinks every row alike, leaving the rows narrower
    room = figure_style.MAX_HEIGHT - TITLE_ROOM * len(heights)
    scale = min(1, room / sum(heights))
    heights = [h * scale for h in heights]
    fig_height = sum(heights) + TITLE_ROOM * len(heights)
    fig = plt.figure(figsize=(figure_style.WIDTH, fig_height))

    axes = {}
    letters = iter(string.ascii_uppercase)
    top = fig_height
    for row, row_images, height in zip(rows, images, heights):
        top -= TITLE_ROOM
        widths = [height * image.shape[1] / image.shape[0] for image in row_images]
        x = (figure_style.WIDTH - sum(widths) - GAP * (len(row) - 1)) / 2
        for name, image, width in zip(row, row_images, widths):
            ax = fig.add_axes([x / figure_style.WIDTH, (top - height) / fig_height,
                               width / figure_style.WIDTH, height / fig_height])
            ax.imshow(image, interpolation='lanczos')
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_color(BORDER_COLOR)
                spine.set_linewidth(0.6)
            ax.set_title(f'{next(letters)}  {TITLES[name]}', loc='left', fontweight='bold', pad=3)
            axes[name] = ax
            x += width + GAP
        top -= height

    for start, end in arrows:
        fig.add_artist(ConnectionPatch(
            xyA=(1, 0.5), coordsA=axes[start].transAxes, xyB=(0, 0.5), coordsB=axes[end].transAxes,
            arrowstyle='-|>', color='#2f3a46', linewidth=1, shrinkA=1, shrinkB=1))

    for extension in ('svg', 'png'):
        fig.savefig(f'{output_stem}.{extension}', dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--screenshots', default='paper/figures/web')
    parser.add_argument('--figures', default='paper/figures', help='output directory')
    args = parser.parse_args()

    for name, rows in FIGURES.items():
        output_stem = os.path.join(args.figures, name)
        draw(rows, ARROWS.get(name, []), args.screenshots, output_stem)
        print(f"Wrote {output_stem}.svg and .png")
