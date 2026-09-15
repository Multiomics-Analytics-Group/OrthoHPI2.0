'''
Venn of the human proteins the parasite-specific TISSUES filter keeps against the ones
HPA single-cell expression (nTPM > cutoff in a cell type of an infected tissue) keeps,
pooled over the parasites as (parasite, protein) pairs. Only tissues HPA covers count.

The pool is each parasite's predicted human interactors (--universe predictions; these
all passed TISSUES already) or the human STRING proteome (--universe proteome, with
--after-deeploc for the part the DeepLoc filter lets through).

Usage: .venv/bin/python scripts/venn_tissues_vs_hpa.py [--universe predictions|proteome]
           [--after-deeploc] [--ntpm 1.0] [--output FILE]
'''
import argparse
import math
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from pipeline import cell_type_annotations, filters, main as pipeline_main
import compare_tissue_filters

HUMAN = 9606
# the app's HPA cell-type cutoff (web_utils.CELL_TYPE_CUTOFFS)
NTPM_CUTOFF = 1.0

# dataviz palette
TISSUES_COLOR = '#2a78d6'
HPA_COLOR = '#eb6834'
BOTH_COLOR = '#4a3aa7'
NEITHER_COLOR = '#c9c8c2'
INK = '#0b0b0b'
INK_SECONDARY = '#52514e'


def hpa_tissue_proteins(config_file, ntpm_cutoff):
    '''{config tissue label: human proteins expressed above the cutoff in some cell type}.'''
    data = cell_type_annotations.read_hpa(config_file=config_file)
    data = cell_type_annotations.map_hpa_data(config_file=config_file, hpa_data=data)
    data = data.dropna(subset=['Gene'])
    data = data[data['nTPM'] > ntpm_cutoff]
    return data.groupby('Tissue')['Gene'].apply(set).to_dict()


def tissues_tissue_proteins(data_dir):
    '''{config tissue label: human proteins the TISSUES filter annotated to it}.'''
    tissues = utils.read_parquet_file(input_file=os.path.join(data_dir, 'tissues_cell_types.parquet'))
    tissues = tissues[tissues['Gene'].str.startswith(f'{HUMAN}.')]
    return tissues.drop_duplicates(['Gene', 'Tissue']).groupby('Tissue')['Gene'].apply(set).to_dict()


def raw_tissues_tissue_proteins(config_file, pool, taxid=HUMAN):
    '''{config tissue label: pool proteins TISSUES annotates to it}, from the raw file.'''
    hosts = utils.read_config(filepath=config_file, field='hosts')
    mapping = utils.read_config(filepath=config_file, field='tissues')
    filename = utils.download_file(url=hosts[taxid]['tissues_url'], data_dir='data/downloads')
    annotations, _ = filters.get_tissues(
        config_file, filename, dict(pool), filters.tissue_cutoff(hosts[taxid], pipeline_main.TISSUE_CUTOFF),
        mapping, taxid)
    by_tissue = {}
    for protein, tissues in annotations.items():
        for tissue in tissues:
            by_tissue.setdefault(tissue.lower(), set()).add(protein)
    return by_tissue


def prediction_universes(data_dir):
    '''{parasite taxid: its predicted human interactors}.'''
    predictions = utils.read_parquet_file(input_file=os.path.join(data_dir, 'predictions.parquet'))
    predictions = predictions[predictions['taxid2'].astype(int) == HUMAN]
    return predictions.groupby(predictions['taxid1'].astype(int))['target'].apply(set).to_dict()


def proteome_universe(config_file, after_deeploc):
    '''The human STRING proteome, optionally narrowed to what the DeepLoc filter keeps.'''
    pool = compare_tissue_filters.host_pools(config_file)
    if after_deeploc:
        return compare_tissue_filters.deeploc_pass(config_file, pool)[HUMAN], pool[HUMAN]
    return set(pool[HUMAN]), pool[HUMAN]


def compare(config_file, data_dir, ntpm_cutoff, universe='predictions', after_deeploc=False):
    '''One row per parasite with the four Venn regions, plus the tissues HPA lacks.'''
    config = utils.read_config(config_file)
    labels = {bto: name.lower() for bto, name in config['tissues'].items()}
    hpa = hpa_tissue_proteins(config_file, ntpm_cutoff)

    if universe == 'proteome':
        proteome, pool = proteome_universe(config_file, after_deeploc)
        tis = raw_tissues_tissue_proteins(config_file, pool)
        targets = {int(taxid): proteome for taxid, parasite in config['parasites'].items()
                   if parasite.get('hosts') is None or HUMAN in parasite['hosts']}
    else:
        tis = tissues_tissue_proteins(data_dir)
        targets = prediction_universes(data_dir)

    rows, skipped = [], {}
    for taxid, parasite in config['parasites'].items():
        if int(taxid) not in targets:
            continue
        tissues = [labels[t] for t in parasite['tissues']]
        covered = [t for t in tissues if t in hpa]
        uncovered = [t for t in tissues if t not in hpa]
        if uncovered:
            skipped[parasite['label']] = uncovered
        if not covered:
            continue

        universe = targets[int(taxid)]
        keep_tis = universe & set().union(*(tis.get(t, set()) for t in covered))
        keep_hpa = universe & set().union(*(hpa[t] for t in covered))
        rows.append({
            'parasite': parasite['label'],
            'taxid': int(taxid),
            'tissues': ', '.join(covered),
            'interactors': len(universe),
            'both': len(keep_tis & keep_hpa),
            'TISSUES only': len(keep_tis - keep_hpa),
            'HPA only': len(keep_hpa - keep_tis),
            'neither': len(universe - keep_tis - keep_hpa),
            '_hpa_only': keep_hpa - keep_tis,
            '_tissues_only': keep_tis - keep_hpa,
            '_tis': keep_tis,
            '_hpa': keep_hpa,
            '_universe': universe,
        })

    return pd.DataFrame(rows), skipped


def pool_pairs(table):
    '''The four regions summed over the parasites, in (parasite, protein) pairs.'''
    totals = table[['interactors', 'both', 'TISSUES only', 'HPA only', 'neither']].sum().to_dict()
    totals['unit'] = '(parasite, human protein) pairs'
    return totals


def pool_proteins(table):
    '''The four regions over unique human proteins, a filter keeping a protein for any parasite.'''
    tis = set().union(*table['_tis'])
    hpa = set().union(*table['_hpa'])
    universe = set().union(*table['_universe'])
    return {'interactors': len(universe), 'both': len(tis & hpa), 'TISSUES only': len(tis - hpa),
            'HPA only': len(hpa - tis), 'neither': len(universe - tis - hpa),
            'unit': 'human proteins, kept for at least one parasite'}


def lens_area(r1, r2, d):
    '''Area shared by two circles of radii r1, r2 whose centres are d apart.'''
    if d >= r1 + r2:
        return 0.0
    if d <= abs(r1 - r2):
        return math.pi * min(r1, r2) ** 2
    a1 = r1 ** 2 * math.acos((d ** 2 + r1 ** 2 - r2 ** 2) / (2 * d * r1))
    a2 = r2 ** 2 * math.acos((d ** 2 + r2 ** 2 - r1 ** 2) / (2 * d * r2))
    k = math.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
    return a1 + a2 - k / 2


def circle_distance(r1, r2, shared):
    '''Centre distance at which the two circles overlap by exactly `shared`.'''
    lo, hi = abs(r1 - r2), r1 + r2
    for _ in range(100):
        mid = (lo + hi) / 2
        if lens_area(r1, r2, mid) > shared:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def draw_venn(totals, ntpm_cutoff, output, universe_label):
    '''Area-proportional Venn: the two filters inside a circle that is the whole pool.'''
    left, both, right = totals['TISSUES only'], totals['both'], totals['HPA only']
    r1 = math.sqrt((left + both) / math.pi)
    r2 = math.sqrt((right + both) / math.pi)
    d = circle_distance(r1, r2, both)
    x_min, x_max = min(-r1, d - r2), max(r1, d + r2)
    # the pool circle, centred on the two filters
    R = math.sqrt(totals['interactors'] / math.pi)
    cx = (x_min + x_max) / 2
    scale = R

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_xlim(cx - 1.6 * R, cx + 1.6 * R)
    ax.set_ylim(-1.55 * R, 1.35 * R)

    ax.add_patch(Circle((cx, 0), R, facecolor=NEITHER_COLOR, alpha=0.45, edgecolor=INK_SECONDARY,
                        linewidth=1.2))
    for x, r, color in ((0, r1, TISSUES_COLOR), (d, r2, HPA_COLOR)):
        ax.add_patch(Circle((x, 0), r, facecolor=color, alpha=0.5, edgecolor=color, linewidth=2))

    def leader(label, xy, xytext):
        ax.annotate(label, xy=xy, xytext=xytext, ha='center', va='center', fontsize=15,
                    color=INK, fontweight='bold',
                    arrowprops={'arrowstyle': '-', 'color': INK_SECONDARY, 'linewidth': 1, 'shrinkB': 2})

    # horizontal extent of each region along the centre line
    regions = [((-r1, d - r2), left, -1), ((d - r2, r1), both, 0), ((r1, d + r2), right, 1)]
    for (a, b), n, side in regions:
        centre, width = (a + b) / 2, b - a
        if width >= 0.25 * scale:
            ax.text(centre, 0, f'{n:,}', ha='center', va='center', fontsize=15, color=INK,
                    fontweight='bold')
        else:
            # too thin for its count: label outside the pool
            leader(f'{n:,}', (centre, 0), (cx + side * 1.3 * R, -0.85 * R))

    # the pool's own count goes in the ring above the filters, or outside if the ring is thin
    inner_top = max(math.sqrt(max(r1 ** 2 - cx ** 2, 0)), math.sqrt(max(r2 ** 2 - (cx - d) ** 2, 0)))
    if R - inner_top >= 0.12 * scale:
        ax.text(cx, (R + inner_top) / 2, f"{totals['neither']:,}", ha='center', va='center',
                fontsize=15, color=INK, fontweight='bold')
    else:
        leader(f"{totals['neither']:,}", (cx, (R + inner_top) / 2), (cx, 1.18 * R))

    handles = [Patch(facecolor=NEITHER_COLOR, alpha=0.45, edgecolor=INK_SECONDARY,
                     label=f'neither: whole pool, {universe_label}'),
               Patch(facecolor=TISSUES_COLOR, alpha=0.5, edgecolor=TISSUES_COLOR,
                     label='TISSUES keeps: ≥ 2.5 in an infected tissue'),
               Patch(facecolor=HPA_COLOR, alpha=0.5, edgecolor=HPA_COLOR,
                     label=f'HPA keeps: nTPM > {ntpm_cutoff:g} in a cell type of an infected tissue')]
    ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.02), ncol=1,
              frameon=False, fontsize=10, labelcolor=INK, handlelength=1.4)

    ax.set_title('Which human proteins each parasite-specific tissue filter keeps\n'
                 f"{totals['interactors']:,} {totals['unit']}; "
                 f"{len(totals['parasites'])} parasites, HPA-covered tissues only; areas proportional to counts",
                 fontsize=11, color=INK, pad=10)

    fig.savefig(output, dpi=200, bbox_inches='tight', facecolor='white')
    return output


def draw_barplot(table, ntpm_cutoff, output, universe_label):
    '''One 100%-stacked bar per parasite over the four regions, sorted by HPA only.'''
    table = table.sort_values('HPA only', ascending=True)
    segments = [('both', 'kept by both', BOTH_COLOR), ('HPA only', 'HPA only (TISSUES drops)', HPA_COLOR),
                ('TISSUES only', 'TISSUES only (HPA drops)', TISSUES_COLOR), ('neither', 'neither', NEITHER_COLOR)]
    share = table[[c for c, _, _ in segments]].div(table['interactors'], axis=0) * 100

    fig, ax = plt.subplots(figsize=(9, 0.28 * len(table) + 1.8))
    left = pd.Series(0.0, index=table.index)
    for column, label, color in segments:
        ax.barh(table['parasite'], share[column], left=left, color=color, label=label,
                height=0.72, edgecolor='white', linewidth=0.6)
        left = left + share[column]
    for y, (n, pct) in enumerate(zip(table['HPA only'], share['HPA only'])):
        ax.text(share['both'].iloc[y] + pct / 2, y, f'{n:,}', ha='center', va='center',
                fontsize=7.5, color='white' if pct > 8 else INK)

    ax.set_xlim(0, 100)
    ax.set_xlabel(f'% of {universe_label}', color=INK_SECONDARY)
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(decimals=0))
    ax.tick_params(axis='y', labelsize=8.5, length=0)
    ax.tick_params(axis='x', colors=INK_SECONDARY)
    for spine in ('top', 'right', 'left'):
        ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_color('#c9c8c2')
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=4, frameon=False, fontsize=8.5)
    ax.set_title(f'Human proteins each parasite-specific tissue filter keeps '
                 f'(TISSUES ≥ 2.5, HPA nTPM > {ntpm_cutoff:g}; HPA-covered tissues only)',
                 fontsize=10, color=INK, pad=28)

    fig.savefig(output, dpi=200, bbox_inches='tight', facecolor='white')
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--universe', choices=['predictions', 'proteome'], default='predictions',
                        help='the human proteins each parasite is judged on (default: predictions)')
    parser.add_argument('--after-deeploc', action='store_true',
                        help='with --universe proteome: only the proteins the DeepLoc filter keeps')
    parser.add_argument('--ntpm', type=float, default=NTPM_CUTOFF, help='HPA nTPM cutoff (default: the app\'s cell-type cutoff)')
    parser.add_argument('--barplot', action='store_true', help='one stacked bar per parasite instead of the Venn')
    parser.add_argument('--proteins', action='store_true',
                        help='count unique human proteins (kept for any parasite) instead of (parasite, protein) pairs')
    parser.add_argument('--output', help='PNG to write (default: snapshots/tissues_vs_hpa_<venn|bars>_<universe>.png)')
    args = parser.parse_args()
    if args.after_deeploc and args.universe != 'proteome':
        parser.error('--after-deeploc only applies to --universe proteome '
                     '(predictions.parquet is already DeepLoc-filtered)')

    universe_label = {'predictions': 'predictions.parquet',
                      'proteome': 'the DeepLoc-passing human proteome' if args.after_deeploc
                      else 'the human STRING proteome'}[args.universe]
    output = args.output or os.path.join(
        'snapshots', f"tissues_vs_hpa_{'bars' if args.barplot else 'venn'}_{args.universe}"
                     f"{'_deeploc' if args.after_deeploc else ''}{'_proteins' if args.proteins else ''}.png")

    table, skipped = compare(args.config, args.data_dir, args.ntpm, args.universe, args.after_deeploc)
    regions = ['interactors', 'both', 'TISSUES only', 'HPA only', 'neither']
    totals = pool_proteins(table) if args.proteins else pool_pairs(table)
    totals['parasites'] = table['parasite'].tolist()

    pd.set_option('display.width', 200)
    print(table[['parasite', 'tissues'] + regions].to_string(index=False))
    print()
    print(f"pooled {totals['unit']}:", {k: int(totals[k]) for k in regions})
    hpa_only = set().union(*table['_hpa_only'])
    tissues_only = set().union(*table['_tissues_only'])
    print(f"unique human proteins HPA keeps and TISSUES drops: {len(hpa_only)}")
    print(f"unique human proteins TISSUES keeps and HPA drops: {len(tissues_only)}")
    if skipped:
        print('\ntissues without HPA data, left out of the comparison:')
        for parasite, tissues in skipped.items():
            print(f"  {parasite}: {', '.join(tissues)}")

    os.makedirs(os.path.dirname(output) or '.', exist_ok=True)
    draw = draw_barplot if args.barplot else draw_venn
    data = table if args.barplot else totals
    print('\nwrote', draw(data, args.ntpm, output, universe_label))


if __name__ == '__main__':
    main()
