'''
Count how many proteins survive each stage of the pipeline, per host-parasite pair, and
draw it as the two branches of the workflow figure. On the host side: proteome ->
expressed in a tissue this parasite infects (TISSUES) -> put by DeepLoc where its niche
reaches -> in a predicted interaction with this parasite. On the
parasite side: proteome -> secreted or on the surface (DeepLoc) -> in an EggNOG group
-> in a predicted interaction with this host. Writes paper/tables/funnel.csv (the counts)
and paper/figures/funnel.pdf/.svg, a butterfly chart: one row per pair, the host stages
nested leftwards and the parasite stages rightwards of the parasite's name, on log axes
unless --linear, which gives each side a linear axis up to its own largest proteome.
'''
import argparse
import copy
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from pipeline import filters, homology, main
from scripts import figure_style

# the orthologous group is left off the host side: every host protein that reaches this
# far is already in one, so the band is zero wide in every row and only costs contrast
HOST_STAGES = ['proteome', 'tissue', 'localization', 'interactor']
PARASITE_STAGES = ['proteome', 'localization', 'eggnog', 'interactor']
# the stage names follow the boxes of the workflow figure
STAGE_LABELS = {
    'proteome': 'STRING proteome',
    'tissue': 'TISSUES filter: expressed in a tissue the parasite infects (host side)',
    'localization': 'DeepLoc filter: in a relevant subcellular location',
    'eggnog': 'in an EggNOG orthologous group (parasite side)',
    'interactor': 'in a predicted host–parasite PPI',
}
# one hue, later stage darker; keyed by stage so both sides agree with the legend. Each
# side draws four of the five, and the two subsets differ, so neighbouring bands within
# a bar are two steps of the ramp apart rather than one
STAGE_COLORS = {
    'proteome': '#e4ecf7',
    'tissue': '#b9d0ec',
    'localization': '#87b0de',
    'eggnog': '#4a7fc1',
    'interactor': '#10315e',
}
# the host-side DeepLoc filter turns on the niche, so the rows are ordered by it; the
# wider-reaching niche first, as in the workflow figure
NICHE_ORDER = ['intracellular', 'extracellular']
NICHE_COLORS = {'intracellular': '#555555', 'extracellular': '#bdbdbd'}
# the colour strips left of a block, as (column, colours, offset in axes fractions)
STRIPS = [('niche', NICHE_COLORS, -0.10), ('group', None, -0.05)]
STRIP_WIDTH = 0.025
# the supplementary table: the columns of the counts under the names a reader needs,
# in the order the figure draws the rows
SUPPLEMENTARY_COLUMNS = {
    'host': 'Host',
    'species': 'Parasite',
    'taxid': 'Parasite taxid',
    'group': 'Taxonomic group',
    'niche': 'Niche in the host',
    'host_proteome': 'Host proteins: STRING proteome',
    'host_tissue': 'Host proteins: expressed in an infected tissue',
    'host_localization': 'Host proteins: in a reachable localisation',
    'host_eggnog': 'Host proteins: in an orthologous group',
    'host_interactor': 'Host proteins: in a predicted PPI',
    'parasite_proteome': 'Parasite proteins: STRING proteome',
    'parasite_localization': 'Parasite proteins: secreted or on the surface',
    'parasite_eggnog': 'Parasite proteins: in an orthologous group',
    'parasite_interactor': 'Parasite proteins: in a predicted PPI',
    'edges': 'Predicted interactions',
}
# a host with fewer pairs still gets a block this many rows tall
MIN_BLOCK = 3
# the layout in inches: margins, the gap between the two columns, and the gap between a
# row's host and parasite bars that holds the parasite's name
LEFT, RIGHT, TOP, BOTTOM = 0.15, 0.1, 1.45, 1.45
COLUMN_GAP, NAME_GAP, BLOCK_GAP = 0.35, 1.05, 4


def load(config_file, data_dir):
    '''Everything the counts read: the proteomes, per host protein its infected-tissue
    labels and the niches reaching it, the eligible pool, which proteins sit in an
    EggNOG group, and the predictions.'''
    downloads = os.path.join(data_dir, 'downloads')
    print("Getting proteins...")
    proteins = main.get_proteins(config_file)
    print("Applying the tissue filter...")
    tissues = filters.apply_tissue_filter(config_file=config_file,
                                          valid_proteins=copy.deepcopy(proteins),
                                          cutoff=main.TISSUE_CUTOFF)
    print("Applying the DeepLoc filter...")
    reachable = filters.apply_deeploc_filter(
        config_file=config_file, valid_proteins=copy.deepcopy(proteins),
        deeploc_dir=os.path.join(data_dir, main.DEEPLOC_ACCURATE_DIR),
        niche_cutoffs=main.DEEPLOC_NICHE_CUTOFFS, default_niche=main.DEEPLOC_DEFAULT_NICHE)

    # eligible_proteins.parquet holds what the secretome, tissue and DeepLoc filters left
    eligible = pd.read_parquet(os.path.join(data_dir, 'eligible_proteins.parquet'))
    eligible['taxid'] = eligible['taxid'].astype(int)
    print("Scanning EggNOG groups...")
    groups = homology.get_eggnog_groups(os.path.join(downloads, '2759_members.tsv.gz'),
                                        eligible['protein'])
    in_group = set().union(*groups.values()) if groups else set()
    predictions = pd.read_parquet(os.path.join(data_dir, 'predictions.parquet'))

    return proteins, tissues, reachable, eligible, in_group, predictions


def count_pairs(config, data):
    '''One row per host-parasite pair: the host proteins left after every host stage,
    the tissue and localization filters being the parasite's own, and the parasite
    proteins left after every parasite stage, the last being those in a predicted
    interaction with this host.'''
    proteins, tissues, reachable, eligible, in_group, predictions = data
    group_order = list(config['parasite_groups'])
    hosts = {t: h['label'].split('(')[1].rstrip(')') for t, h in config['hosts'].items()}
    # niche first, so that the host side's two regimes sit together; congeners share a
    # niche, so ordering by group and name within it keeps them adjacent as before
    def order(item):
        return (NICHE_ORDER.index(item[1]['niche']), group_order.index(item[1]['group']),
                item[1]['label'])

    rows = []
    for host, host_name in hosts.items():
        prefix = f'{host}.'
        for taxid, parasite in sorted(config['parasites'].items(), key=order):
            if host not in parasite['hosts']:
                continue
            labels = {config['tissues'][t] for t in parasite['tissues']}
            expressed = {p for p, ts in tissues.items() if p.startswith(prefix) and labels & set(ts)}
            reach = expressed & reachable[parasite['niche']]
            pool = set(eligible.loc[eligible['taxid'] == taxid, 'protein'])
            edges = predictions[(predictions['taxid1'] == str(taxid)) & (predictions['taxid2'] == str(host))]
            rows.append({'host': host_name, 'taxid': taxid, 'species': parasite['label'],
                         'group': parasite['group'], 'niche': parasite['niche'],
                         'host_proteome': len(proteins[host]), 'host_tissue': len(expressed),
                         'host_localization': len(reach), 'host_eggnog': len(reach & in_group),
                         'host_interactor': edges['target'].nunique(),
                         'parasite_proteome': len(proteins[taxid]), 'parasite_localization': len(pool),
                         'parasite_eggnog': len(pool & in_group),
                         'parasite_interactor': edges['source'].nunique(),
                         'edges': len(edges)})

    return pd.DataFrame(rows)


def draw_side(ax, rows, side, stages, mirrored, limit=None):
    '''Nested bars of the proteins each stage keeps on one side, one row per pair, on a
    log axis that grows leftwards when `mirrored`, or a linear one up to `limit`.'''
    y = range(len(rows))
    for stage in stages:
        # the bars are nested, so each one's outline draws the boundary of its band
        # onto the lighter band under it, visible even where the band is a point wide
        if limit is None:
            # the log axis starts at 1, so a bar starts there too; a stage with nothing
            # left draws no bar
            ax.barh(y, rows[f'{side}_{stage}'].astype(float).clip(lower=1) - 1, left=1,
                    color=STAGE_COLORS[stage], height=0.72, edgecolor='white', linewidth=0.5)
        else:
            ax.barh(y, rows[f'{side}_{stage}'].astype(float), color=STAGE_COLORS[stage],
                    height=0.72, edgecolor='white', linewidth=0.5)
    ax.set_yticks(list(y))
    ax.set_ylim(max(len(rows), MIN_BLOCK) - 0.4, -0.6)
    if limit is None:
        ax.set_xscale('log')
        ax.set_xlim((1e5, 1) if mirrored else (1, 1e5))
        ax.set_xticks([1, 10, 100, 1000, 10000, 100000])
        ax.set_xticklabels(['1', '', '100', '', '10,000', ''])
        ax.minorticks_off()
    else:
        ax.set_xlim((limit, 0) if mirrored else (0, limit))
        ax.set_xticks([0, limit // 2, limit])
        # the midpoint keeps its gridline but not a label, which the narrow column of a
        # block has no room for
        ax.set_xticklabels(['0', '', f'{limit:,}'])
        # the outermost label stands at the edge of the figure, so it reads inwards
        ax.get_xticklabels()[-1].set_ha('left' if mirrored else 'right')
    ax.tick_params(axis='y', length=0)
    ax.grid(axis='x', color='#e6e6e6', linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ('top', 'right', 'left'):
        ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_color('#999999')


def add_strips(ax, rows, group_colors):
    '''Left of the bars, one colour strip per kind: the niche outermost, the taxonomic
    group beside it. A strip spans a run of neighbouring rows sharing a value, so a
    group split across the two niches gets a strip in each.'''
    for column, colors, offset in STRIPS:
        colors = colors or group_colors
        values = list(rows[column])
        start = 0
        for i in range(1, len(values) + 1):
            if i == len(values) or values[i] != values[start]:
                ax.add_patch(Rectangle((offset, start - 0.4), STRIP_WIDTH, i - start - 0.2,
                                       color=colors[values[start]], linewidth=0,
                                       transform=ax.get_yaxis_transform(), clip_on=False))
                start = i


def draw_block(fig, rows, group_colors, x, y_top, row, width, limits):
    '''One host's pairs as a butterfly: the host bars growing left and the parasite
    bars growing right of the parasite names, which sit in the gap between them, and
    the niche and group strips down the far left. Returns the pair of axes.'''
    height = max(len(rows), MIN_BLOCK) * row
    bar_width = (width - NAME_GAP) / 2
    figure_width, figure_height = fig.get_size_inches()
    host = fig.add_axes([x / figure_width, (y_top - height) / figure_height,
                         bar_width / figure_width, height / figure_height])
    parasite = fig.add_axes([(x + bar_width + NAME_GAP) / figure_width, (y_top - height) / figure_height,
                             bar_width / figure_width, height / figure_height])
    draw_side(host, rows, 'host', HOST_STAGES, mirrored=True, limit=limits.get('host'))
    draw_side(parasite, rows, 'parasite', PARASITE_STAGES, mirrored=False,
              limit=limits.get('parasite'))
    host.set_yticklabels([])
    # the name is centred in the gap, so its anchor is pushed half the gap from the axis
    parasite.set_yticklabels([figure_style.short_name(s) for s in rows['species']],
                             style='italic', ha='center')
    parasite.tick_params(axis='y', pad=NAME_GAP / 2 * 72)
    add_strips(host, rows, group_colors)
    host.set_title(rows['host'].iloc[0], loc='left', fontweight='bold')

    return host, parasite


def draw(table, group_colors, output_stem, linear=False):
    '''The host with the most pairs in the left column, the other hosts stacked in the
    right, the stage legend above and the group legend below.'''
    figure_style.apply()
    # a linear axis runs to the largest proteome of its own side, the log one to 1e5 on
    # both, so the two sides keep a common scale only on the log axis
    limits = ({side: int(table[f'{side}_proteome'].max()) for side in ('host', 'parasite')}
              if linear else {})
    hosts = list(dict.fromkeys(table['host']))
    sizes = {h: max((table['host'] == h).sum(), MIN_BLOCK) for h in hosts}
    first = max(hosts, key=sizes.get)
    others = [h for h in hosts if h != first]
    n_rows = max(sizes[first], sum(sizes[h] for h in others) + BLOCK_GAP * (len(others) - 1))
    # the rows shrink only when a page would not hold them at the shared row height
    row = min(figure_style.ROW, (figure_style.MAX_HEIGHT - TOP - BOTTOM) / n_rows)
    fig = plt.figure(figsize=(figure_style.WIDTH, row * n_rows + TOP + BOTTOM))
    width = (figure_style.WIDTH - LEFT - RIGHT - COLUMN_GAP) / 2
    y_top = fig.get_figheight() - TOP

    blocks = [(first, LEFT, y_top)]
    y = y_top
    for host in others:
        blocks.append((host, LEFT + width + COLUMN_GAP, y))
        y -= (sizes[host] + BLOCK_GAP) * row
    axes = {}
    for host, x, y in blocks:
        rows = table[table['host'] == host].reset_index(drop=True)
        axes[host] = draw_block(fig, rows, group_colors, x, y, row, width, limits)
    for host in (first, others[-1]):
        axes[host][0].set_xlabel('Host proteins')
        axes[host][1].set_xlabel('Parasite proteins')

    handles = [Patch(color=color, label=STAGE_LABELS[stage]) for stage, color in STAGE_COLORS.items()]
    fig.legend(handles=handles, title='Proteins left after each stage', loc='upper left', ncol=1, frameon=False,
               bbox_to_anchor=(LEFT / figure_style.WIDTH, 0.995), handlelength=1.2,
               alignment='left', title_fontsize=figure_style.FONT)
    figure_style.stack_legends(fig, left=LEFT / figure_style.WIDTH, blocks=[
        ('Parasite niche in the host',
         [Patch(color=NICHE_COLORS[niche], label=niche) for niche in NICHE_ORDER]),
        ('Taxonomic group', [Patch(color=color, label=group) for group, color in group_colors.items()])])
    for extension in ('pdf', 'svg'):
        fig.savefig(f'{output_stem}.{extension}')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--table', default='paper/tables/funnel.csv')
    parser.add_argument('--supplementary-table', default='paper/tables/funnel_supplementary.csv',
                        help='the same rows under readable headers')
    parser.add_argument('--figure', default='paper/figures/funnel', help='output path without extension')
    parser.add_argument('--linear', action='store_true', help='linear axes instead of the log ones')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    pairs = count_pairs(config, load(args.config, args.data_dir))
    pairs.to_csv(args.table, index=False)
    supplementary = pairs[list(SUPPLEMENTARY_COLUMNS)].rename(columns=SUPPLEMENTARY_COLUMNS)
    supplementary.to_csv(args.supplementary_table, index=False)
    print(f"Wrote {args.table} and {args.supplementary_table} ({len(pairs)} pairs)")
    draw(pairs, config['parasite_groups'], args.figure, linear=args.linear)
    print(f"Wrote {args.figure} (.pdf and .svg)")
