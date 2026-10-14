'''
Screenshots of the running app for the web interface figures, written to
paper/figures/web/ for scripts/build_web_interface_figure.py.

  home              the home page's hosts and parasites
  multi_host        the multi-host parasites page for C. muris
  parasite_overlap  the parasites-of-a-host page for human, the parasites ordered by niche
                    as in the paper's similarity figure, the P. falciparum - P. vivax
                    cell clicked
  network, alphafold, cell_types, enrichment
                    the network page: P. falciparum in human at confidence 0.7, liver and
                    20 nTPM, a link clicked for its AlphaFold models

Needs the app running and Playwright with Chromium, both in the streamlit environment:

  streamlit run app/OrthoHPI_Home.py --server.headless true --server.port 8599 --theme.base light
  ~/.venvs/streamlit/bin/python scripts/capture_web_screenshots.py
'''
import argparse
import os
import re

import numpy as np
from PIL import Image
from playwright.sync_api import sync_playwright

VIEWPORT = {'width': 1500, 'height': 9000}
# the centred column of the network page, in CSS pixels
COLUMN = (390, 1110)
HIDE_CAPTIONS = '[data-testid=stCaptionContainer]{display:none !important}'
# the network component keeps its vis.js network inside a closure; the served page is
# patched to hand it out, so a link can be found and clicked
NETWORK_PAGE = re.compile(r'/component/ppi_network[^/]*/index\.html')
# the link whose midpoint lies furthest from any node, so the click lands on it alone
FIND_LINK = '''() => {
  const net = window.__getNet();
  const pos = net.getPositions();
  const box = document.getElementById('mynetwork').getBoundingClientRect();
  const dom = id => net.canvasToDOM(pos[id]);
  let best = null;
  for (const e of net.body.data.edges.get()) {
    const a = dom(e.from), b = dom(e.to);
    const m = {x: (a.x + b.x) / 2, y: (a.y + b.y) / 2};
    let near = Math.hypot(a.x - b.x, a.y - b.y) / 2;
    for (const id in pos) { const p = dom(id); near = Math.min(near, Math.hypot(p.x - m.x, p.y - m.y)); }
    if (!best || near > best.near) best = {near, x: m.x + box.left, y: m.y + box.top};
  }
  return best;
}'''
# the centre of the heatmap cell of two parasites
FIND_CELL = '''([a, b]) => {
  const gd = document.querySelector('.js-plotly-plot');
  const L = gd._fullLayout;
  const names = L.yaxis.ticktext || gd.data.find(t => t.type === 'heatmap').y;
  const i = names.findIndex(n => n.includes(a)), j = names.findIndex(n => n.includes(b));
  const r = gd.getBoundingClientRect();
  return {x: r.left + L.xaxis._offset + L.xaxis.d2p(j), y: r.top + L.yaxis._offset + L.yaxis.d2p(i)};
}'''


def expose_network(route):
    response = route.fetch()
    body = response.text().replace(
        'var network = null;', 'var network = null; window.__getNet = function () { return network; };', 1)
    route.fulfill(response=response, body=body)


def wait(page):
    '''Until Streamlit has finished the run.'''
    page.wait_for_timeout(1000)
    page.wait_for_function("() => !document.querySelector('[data-testid=stStatusWidget]')", timeout=240000)
    page.wait_for_timeout(1000)


def choose(page, label, text):
    page.get_by_label(label).first.click()
    page.keyboard.type(text)
    page.wait_for_timeout(300)
    page.keyboard.press('Enter')
    wait(page)


def slide(page, index, steps):
    page.get_by_role('slider').nth(index).focus()
    for _ in range(steps):
        page.keyboard.press('ArrowRight')
    wait(page)


def hide_captions(page):
    page.add_style_tag(content=HIDE_CAPTIONS)
    page.wait_for_timeout(1500)


def column_clip(top, bottom, pad=12):
    '''The network page's column from the top of one element to the bottom of another.'''
    t, b = top.bounding_box(), bottom.bounding_box()
    return {'x': COLUMN[0], 'y': t['y'] - pad, 'width': COLUMN[1] - COLUMN[0],
            'height': b['y'] + b['height'] + 2 * pad - t['y']}


def collapse(image, keep=60, margin=60, tol=6):
    '''Shrinks each band of blank rows taller than `keep` to `keep`, testing all but the
    `margin` at either side.'''
    pixels = np.asarray(image.convert('RGB')).astype(int)
    core = pixels[:, margin:pixels.shape[1] - margin]
    blank = (core.max(axis=1) - core.min(axis=1)).max(axis=1) <= tol
    rows, run = [], []
    for i, is_blank in enumerate(blank):
        if is_blank:
            run.append(i)
            continue
        rows += run[:keep // 2] + run[-(keep // 2):] if len(run) > keep else run
        rows.append(i)
        run = []
    rows += run[:keep]

    return Image.fromarray(pixels[rows].astype(np.uint8))


def side_by_side(left, right, gap=40, top=0):
    out = Image.new('RGB', (left.width + gap + right.width, max(left.height, top + right.height)), 'white')
    out.paste(left, (0, 0))
    out.paste(right, (left.width + gap, top))

    return out


def network_page(page, base, out):
    page.goto(f'{base}/Predicted_Host-Parasite_PPIs')
    wait(page)
    choose(page, 'Select a host', 'Human')
    choose(page, 'Select a parasite', 'falciparum')
    slide(page, 0, 35)   # confidence 0.35 -> 0.70
    page.get_by_label('Select tissues to filter').first.click()
    page.keyboard.type('liver')
    page.wait_for_timeout(300)
    page.keyboard.press('Enter')
    page.keyboard.press('Escape')
    wait(page)
    slide(page, 1, 19)   # cell-type expression 1 -> 20 nTPM

    # the AlphaFold dialog, before the captions go: the dialog's line on the STRING link is one
    frame = next(f for f in page.frames if 'ppi_network' in f.url)
    page.wait_for_timeout(4000)  # the network fits itself once the physics settle
    link = frame.evaluate(FIND_LINK)
    network = page.locator('iframe[src*="ppi_network"]')
    box = network.bounding_box()
    page.mouse.click(box['x'] + link['x'], box['y'] + link['y'])
    dialog = page.get_by_role('dialog')
    dialog.wait_for(timeout=120000)
    wait(page)
    page.wait_for_timeout(5000)  # 3Dmol draws after the run ends
    dialog.screenshot(path=f'{out}/alphafold.png')
    collapse(Image.open(f'{out}/alphafold.png')).save(f'{out}/alphafold.png')
    page.keyboard.press('Escape')
    wait(page)

    hide_captions(page)
    page.mouse.move(box['x'] + link['x'], box['y'] + link['y'])
    page.wait_for_timeout(1000)
    page.screenshot(path=f'{out}/network_controls.png', clip=column_clip(
        page.get_by_text('Predicted interactions per infected organ'),
        page.get_by_label('Select cell types to filter').first))
    page.screenshot(path=f'{out}/network_graph.png', clip=column_clip(
        page.get_by_text('Extracellular', exact=True).first, network))
    side_by_side(Image.open(f'{out}/network_controls.png'),
                 collapse(Image.open(f'{out}/network_graph.png'), margin=40)).save(f'{out}/network.png')
    for piece in ('network_controls', 'network_graph'):
        os.remove(f'{out}/{piece}.png')

    page.screenshot(path=f'{out}/cell_types.png', clip=column_clip(
        page.get_by_role('heading', name='Cell types expressing the host proteins'),
        page.locator('[data-testid=stPlotlyChart]').first))
    heading = page.get_by_role('heading', name='Enriched biological processes', exact=True).bounding_box()
    top = heading['y'] - 12
    # down to the end of the first chart under the heading, the ranked processes
    chart = next(box for box in (c.bounding_box() for c in page.locator('[data-testid=stPlotlyChart]').all())
                 if box and box['y'] > heading['y'])
    bottom = chart['y'] + chart['height'] + 12
    page.screenshot(path=f'{out}/enrichment.png', clip={'x': COLUMN[0], 'y': top,
                                               'width': COLUMN[1] - COLUMN[0], 'height': bottom - top})


def home_page(page, base, out):
    page.goto(base)
    wait(page)
    hide_captions(page)
    top = page.get_by_role('heading', name='Hosts and parasites').bounding_box()
    bottom = page.get_by_text('Host:', exact=True).first.bounding_box()
    page.screenshot(path=f'{out}/home.png', clip={'x': 70, 'y': top['y'] - 10, 'width': 1360,
                                                  'height': bottom['y'] + bottom['height'] + 14 - top['y']})


def multi_host_page(page, base, out):
    page.goto(f'{base}/Hosts_of_a_Parasite')
    wait(page)
    choose(page, 'Select a parasite to compare', 'Cryptosporidium muris')
    hide_captions(page)
    top = page.get_by_role('heading', name='Predicted interactions per infected organ').bounding_box()
    bottom = page.get_by_text('missing links, one row each').bounding_box()
    page.screenshot(path=f'{out}/multi_host.png', clip={'x': 70, 'y': top['y'] - 10, 'width': 1360,
                                               'height': bottom['y'] + bottom['height'] + 14 - top['y']})


def parasites_page(page, base, out):
    page.goto(f'{base}/Compare_Parasites')
    wait(page)
    choose(page, 'Select a host', 'Human')
    page.get_by_role('radiogroup').first.get_by_text('intracellular / extracellular').click()
    wait(page)
    hide_captions(page)
    heading = page.get_by_role('heading', name='Overlap of the host interactors').bounding_box()
    chart = page.locator('.js-plotly-plot').first.bounding_box()
    page.screenshot(path=f'{out}/overlap_heatmap.png', clip={
        'x': chart['x'] - 10, 'y': heading['y'] - 10, 'width': max(chart['width'], heading['width']) + 20,
        'height': chart['y'] + chart['height'] - heading['y'] + 20})
    cell = page.evaluate(FIND_CELL, ['falciparum', 'vivax'])
    page.mouse.click(cell['x'], cell['y'])
    dialog = page.get_by_role('dialog')
    dialog.wait_for(timeout=120000)
    wait(page)
    dialog.screenshot(path=f'{out}/overlap_dialog.png')

    heatmap, shared = Image.open(f'{out}/overlap_heatmap.png'), Image.open(f'{out}/overlap_dialog.png')
    # the dialog at most of the heatmap's height, level with its grid
    scale = 0.62 * heatmap.height / shared.height
    shared = shared.resize((int(shared.width * scale), int(shared.height * scale)), Image.LANCZOS)
    side_by_side(heatmap, shared, gap=30, top=int(0.12 * heatmap.height)).save(f'{out}/parasite_overlap.png')
    for piece in ('overlap_heatmap', 'overlap_dialog'):
        os.remove(f'{out}/{piece}.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--url', default='http://localhost:8599')
    parser.add_argument('--out', default='paper/figures/web')
    parser.add_argument('--only', nargs='+', help='capture these pages alone, by function name')
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    captures = [home_page, multi_host_page, parasites_page, network_page]
    CAPTURES = [c for c in captures if not args.only or c.__name__ in args.only]
    with sync_playwright() as p:
        browser = p.chromium.launch()
        for capture in CAPTURES:
            page = browser.new_page(viewport=VIEWPORT, device_scale_factor=2)
            page.route(NETWORK_PAGE, expose_network)
            capture(page, args.url, args.out)
            page.close()
            print(f'{capture.__name__}: done')
        browser.close()
    print(f'Wrote the screenshots to {args.out}')
