'''
The GO enrichment panel of the conference poster: the processes enriched among the host
proteins of the network page's T. cruzi example (human, confidence 0.7, heart), read from
the app's ranked dot plot and drawn wide and short to fit under the AlphaFold panel.
Written to poster/figures/enrichment_poster.svg.

Needs the app running and Playwright with Chromium, both in the streamlit environment:

  streamlit run app/OrthoHPI_Home.py --server.headless true --server.port 8599 --theme.base light
  ~/.venvs/streamlit/bin/python poster/capture_poster_enrichment.py
'''
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
from playwright.sync_api import sync_playwright

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
from capture_web_screenshots import VIEWPORT, choose, slide, wait  # noqa: E402

# the traces of the first chart under the enrichment heading: the host proteins' dot plot
READ_CHART = '''() => {
  const heading = [...document.querySelectorAll('h1,h2,h3,h4')]
    .find(h => h.textContent.trim() === 'Enriched biological processes');
  const top = heading.getBoundingClientRect().top;
  const chart = [...document.querySelectorAll('[data-testid=stPlotlyChart] .js-plotly-plot')]
    .find(c => c.getBoundingClientRect().top > top);
  const t = chart.data[0];
  return {terms: t.customdata.map(c => c[0]), odds: t.x, n: t.marker.size,
          significance: t.marker.color,
          scale: chart.layout.coloraxis.colorscale};
}'''
TOP_N = 6
FONT = 'Arial'
TEXT = '#3c434c'
GRID = '#e5e7eb'
# poster millimetres to inches
MM = 1 / 25.4


def read_enrichment(base):
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(viewport=VIEWPORT)
        page.goto(f'{base}/Predicted_Host-Parasite_PPIs')
        wait(page)
        choose(page, 'Select a host', 'Human')
        choose(page, 'Select a parasite', 'cruzi')
        slide(page, 0, 35)   # confidence 0.35 -> 0.70
        # the tissue multiselect, by its input, as its label also names the help icon
        page.locator('[data-testid=stMultiSelect]').nth(0).locator('input').click()
        page.keyboard.type('heart')
        page.wait_for_timeout(500)
        page.keyboard.press('Enter')
        page.keyboard.press('Escape')
        wait(page)
        page.get_by_role('heading', name='Enriched biological processes', exact=True).wait_for()
        page.wait_for_timeout(3000)
        data = page.evaluate(READ_CHART)
        browser.close()
    return data


def values(column):
    '''Plotly hands numeric arrays out typed, with the plain list beside them.'''
    return column['_inputArray'] if isinstance(column, dict) else column


def draw(data, path, width_mm, height_mm):
    data = {k: values(v) if k != 'scale' else v for k, v in data.items()}
    rows = sorted(zip(data['terms'], data['odds'], data['n'], data['significance']),
                  key=lambda r: -r[3])[:TOP_N]
    rows.sort(key=lambda r: r[1])
    terms, odds, n, sig = (list(c) for c in zip(*rows))
    terms = [t[0].upper() + t[1:] for t in terms]
    cmap = LinearSegmentedColormap.from_list('go', [c for _, c in data['scale']])
    every = data['significance']

    plt.rcParams.update({'font.family': FONT, 'font.size': 16, 'text.color': TEXT,
                         'axes.labelcolor': TEXT, 'xtick.color': TEXT, 'ytick.color': TEXT,
                         'svg.fonttype': 'none'})
    fig, ax = plt.subplots(figsize=(width_mm * MM, height_mm * MM))
    points = ax.scatter(odds, range(len(terms)), s=[v * 22 for v in n], c=sig, cmap=cmap,
                        vmin=min(every), vmax=max(every), edgecolors='white', linewidths=1.2,
                        zorder=3)
    ax.set_yticks(range(len(terms)), terms)
    ax.set_xscale('log')
    ax.set_xlabel('Odds ratio')
    ax.xaxis.set_major_locator(LogLocator(subs=(1, 2, 5)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(min(odds) / 1.25, max(odds) * 1.25)
    ax.set_ylim(-0.6, len(terms) - 0.4)
    ax.grid(axis='y', color=GRID, linewidth=1, zorder=0)
    ax.tick_params(axis='y', length=0, pad=8)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color('#9ca3af')
    bar = fig.colorbar(points, ax=ax, fraction=0.025, pad=0.02, aspect=12)
    bar.set_label('−log$_{10}$ FDR')
    bar.outline.set_visible(False)
    fig.tight_layout(pad=0.3)
    fig.savefig(path, transparent=True)
    print(f'Wrote {path}: the {len(terms)} most significant of {len(every)} processes')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--url', default='http://localhost:8599')
    parser.add_argument('--out', default='poster/figures/enrichment_poster.svg')
    parser.add_argument('--width', type=float, default=320, help='panel width in mm')
    parser.add_argument('--height', type=float, default=72, help='panel height in mm')
    parser.add_argument('--cached', action='store_true',
                        help='redraw from the data of the last run instead of the app')
    args = parser.parse_args()

    cache = os.path.splitext(args.out)[0] + '.json'
    if args.cached:
        with open(cache) as f:
            data = json.load(f)
    else:
        data = read_enrichment(args.url)
        with open(cache, 'w') as f:
            json.dump(data, f)
    draw(data, args.out, args.width, args.height)
