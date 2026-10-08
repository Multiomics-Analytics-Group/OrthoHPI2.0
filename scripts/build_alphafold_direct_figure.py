'''
The AlphaFold 3 Server models of the direct-interaction benchmark
(scripts/build_alphafold_direct_jobs.py): for each link whose predicted model and human
reference have both run, the interface confidence of the predicted parasite-host
complex against that of the reference, as ipTM (A) and ipSAE (B), coloured by the
link's confidence score, circles where the reference shares the predicted partner. Above the reference's cutoff and below the predicted model's
are transfers that fail; below both, complexes AlphaFold does not model. Runs on
whatever results are in so far. Writes paper/figures/alphafold_direct.svg/.png and
paper/tables/alphafold_direct.csv.
'''
import argparse
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts import alphafold_results, figure_style
from scripts.build_alphafold_figure import CONFIDENT, FAILED

# ipSAE above which an interface counts as found, as in the benchmark's summaries
IPSAE_FOUND = 0.1
METRICS = [('iptm', 'ipTM', FAILED), ('ipsae', 'ipSAE', IPSAE_FOUND)]
# light to dark blue for low to high confidence score
SCORE_COLORS = LinearSegmentedColormap.from_list('score', ['#86b6ef', '#2a78d6', '#0d366b'])
# references sharing the predicted partner, and those of another human pair of the groups
PARTNER_MARKERS = {True: 'o', False: 'D'}
PARTNER_LABELS = {True: 'Reference shares the predicted partner',
                  False: 'Reference is another human pair of the groups'}


def scores(jobs, results):
    '''Each link with both its models run: the best-ranked model's confidences of each.'''
    rows = []
    for job in jobs.to_dict('records'):
        folder = alphafold_results.job_folder(job['job'], results)
        if folder is not None:
            best = alphafold_results.models(folder)[0]
            rows.append({'link': job['link'], 'kind': job['kind'], 'name A': job['name A'],
                         'name B': job['name B'], 'iptm': best['iptm'], 'ipsae': best['ipsae']})
    models = pd.DataFrame(rows)
    links = jobs[jobs['kind'] == 'predicted'][['link', 'parasite', 'bin', 'score',
                                               'same partner']]
    wide = models.pivot(index='link', columns='kind',
                        values=['iptm', 'ipsae', 'name A', 'name B'])
    wide.columns = [f'{metric} {kind}' for metric, kind in wide.columns]
    wide = links.merge(wide.reset_index(), on='link').dropna(subset=['iptm predicted',
                                                                     'iptm reference'])
    for metric in ('iptm', 'ipsae'):
        for kind in ('predicted', 'reference'):
            wide[f'{metric} {kind}'] = wide[f'{metric} {kind}'].astype(float)

    return wide


def outcome(row, metric, cutoff):
    pred, ref = row[f'{metric} predicted'] >= cutoff, row[f'{metric} reference'] >= cutoff
    return ('both' if pred and ref else 'transfer fails' if ref else
            'only predicted' if pred else 'neither')


def draw(links, out):
    figure_style.apply()
    fig, axes = plt.subplots(1, 2, figsize=(figure_style.WIDTH, 4.0))
    norm = Normalize(0.35, 1.0)
    for ax, (metric, label, cutoff), letter in zip(axes, METRICS, 'AB'):
        ax.axvline(cutoff, color='#b4b2a9', lw=0.6, ls='--', zorder=0)
        ax.axhline(cutoff, color='#b4b2a9', lw=0.6, ls='--', zorder=0)
        ax.plot([0, 1], [0, 1], color='#d3d1c7', lw=0.6, zorder=0)
        for same, rows in links.groupby('same partner'):
            ax.scatter(rows[f'{metric} reference'], rows[f'{metric} predicted'], s=26,
                       marker=PARTNER_MARKERS[same], c=rows['score'], cmap=SCORE_COLORS,
                       norm=norm, edgecolors='white', lw=0.5, zorder=3)
        counts = links.apply(outcome, axis=1, args=(metric, cutoff)).value_counts()
        # each quadrant's count in its corner
        for (x, y, ha, va), name in zip([(0.98, 0.98, 'right', 'top'), (0.98, 0.02, 'right', 'bottom'),
                                         (0.02, 0.02, 'left', 'bottom'), (0.02, 0.98, 'left', 'top')],
                                        ['both', 'transfer fails', 'neither', 'only predicted']):
            ax.text(x, y, f'{name}: {counts.get(name, 0)}', transform=ax.transAxes, ha=ha,
                    va=va, fontsize=figure_style.FONT, color='#5f5e5a')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect('equal')
        ax.set_xlabel(f'{label} of the human reference')
        ax.set_ylabel(f'{label} of the predicted link')
        ax.set_title(letter, loc='left', fontweight='bold')
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
    handles = [Line2D([], [], marker=PARTNER_MARKERS[same], color='#5f5e5a', ls='',
                      markersize=5, label=f"{PARTNER_LABELS[same]} "
                                          f"(n = {(links['same partner'] == same).sum()})")
               for same in (True, False)]
    fig.tight_layout(rect=(0, 0, 0.9, 1))
    bottom = figure_style.stack_legends(fig, [('', handles, 2)], left=0.06)
    fig.subplots_adjust(bottom=(bottom + 0.5) / fig.get_figheight())
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=SCORE_COLORS),
                       cax=fig.add_axes([0.92, 0.25, 0.015, 0.5]))
    bar.set_label('Confidence score of the link')
    bar.outline.set_visible(False)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(f'{out}.svg')
    fig.savefig(f'{out}.png', dpi=300)


def build(jobs_file, results, out, table):
    jobs = pd.read_csv(jobs_file)
    links = scores(jobs, results)
    if links.empty:
        raise SystemExit(f'no link with both its models in {results}')
    draw(links, out)
    for metric, _, cutoff in METRICS:
        links[f'outcome {metric}'] = links.apply(outcome, axis=1, args=(metric, cutoff))
    links.to_csv(table, index=False)
    print(f"{len(links)} of {jobs['link'].nunique()} links with both models run")
    for metric, _, _ in METRICS:
        print(pd.crosstab(links[f'outcome {metric}'], links['same partner'], margins=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--jobs', default='paper/alphafold/direct.csv')
    parser.add_argument('--results', default=alphafold_results.RESULTS,
                        help='where the server downloads are unpacked, as folds_*')
    parser.add_argument('--table', default='paper/tables/alphafold_direct.csv')
    parser.add_argument('--figure', default='paper/figures/alphafold_direct',
                        help='output path without extension')
    args = parser.parse_args()
    build(args.jobs, args.results, args.figure, args.table)
