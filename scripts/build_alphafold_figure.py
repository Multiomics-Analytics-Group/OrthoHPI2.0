'''
The AlphaFold 3 Server models of the case-study links set against their references and
controls (scripts/build_alphafold_jobs.py). Each pair is a row: the interface confidence
(ipTM) of the five models of the predicted parasite-host complex, the host-host
complex it was transferred from, and the parasite protein against random host proteins.
Writes paper/figures/alphafold.svg/.png and paper/tables/alphafold.csv.
'''
import argparse
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts import alphafold_results, figure_style

KINDS = ['reference', 'predicted', 'control']
KIND_LABELS = {'reference': 'Host–host reference', 'predicted': 'Predicted link',
               'control': 'Random host protein'}
KIND_COLORS = {'reference': '#eb6834', 'predicted': '#2a78d6', 'control': '#1baf7a'}
KIND_MARKERS = {'reference': 's', 'predicted': 'o', 'control': '^'}
STUDY_LABELS = {'lb': r'$\it{L.\,braziliensis}$', 'pf': r'$\it{P.\,falciparum}$'}
# ipTM above which a complex is confidently predicted, and below which it likely failed
CONFIDENT, FAILED = 0.8, 0.6


def scores(jobs, results):
    '''The confidences of every model of every job, ranked as the server ranks them.'''
    rows = []
    for job in jobs.to_dict('records'):
        folder = alphafold_results.job_folder(job['job'], results)
        if folder is None:
            raise SystemExit(f"no models of {job['job']} in {results}")
        for model in alphafold_results.models(folder):
            rows.append({'job': job['job'], 'case study': job['case study'],
                         'pair': job['pair'], 'kind': job['kind'], 'partner': job['name B'],
                         **model})

    return pd.DataFrame(rows)


def draw(models, out):
    figure_style.apply()
    pairs = list(dict.fromkeys(models['pair']))
    offsets = {'reference': -0.25, 'predicted': 0, 'control': 0.25}
    fig, ax = plt.subplots(figsize=(figure_style.WIDTH * 0.62,
                                    len(pairs) * 3 * figure_style.ROW + 1.1))
    ax.axvspan(0, FAILED, color='#f0efeb', zorder=0, lw=0)
    ax.axvline(CONFIDENT, color='#b4b2a9', lw=0.6, ls='--', zorder=0)
    for i, pair in enumerate(pairs):
        for kind in KINDS:
            rows = models[(models['pair'] == pair) & (models['kind'] == kind)]
            for job, group in rows.groupby('job', sort=False):
                y = i + offsets[kind]
                rest = group[group['model'] > 0]
                ax.scatter(rest['iptm'], [y] * len(rest), s=14, marker=KIND_MARKERS[kind],
                           facecolors='none', edgecolors=KIND_COLORS[kind], lw=0.7, zorder=2)
                top = group[group['model'] == 0]
                ax.scatter(top['iptm'], [y], s=34, marker=KIND_MARKERS[kind],
                           color=KIND_COLORS[kind], edgecolors='white', lw=0.6, zorder=3)
                if kind == 'control':
                    ax.annotate(group['partner'].iloc[0], (group['iptm'].max(), y),
                                xytext=(5, 0), textcoords='offset points', va='center',
                                fontsize=figure_style.FONT, color='#5f5e5a')
        if i:
            ax.axhline(i - 0.5, color='#d3d1c7', lw=0.5)
    study = dict(zip(models['pair'], models['case study']))
    ax.set_yticks(range(len(pairs)),
                  [f'{p}\n({STUDY_LABELS[study[p]]})' for p in pairs], fontsize=figure_style.FONT)
    ax.set_ylim(len(pairs) - 0.5, -0.5)
    ax.set_xlim(0, 1)
    ax.set_xlabel('Interface confidence (ipTM)')
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis='y', length=0)

    handles = [Line2D([], [], marker=KIND_MARKERS[k], color=KIND_COLORS[k], ls='',
                      markersize=5, label=KIND_LABELS[k]) for k in KINDS]
    handles.append(Line2D([], [], marker='o', color='#5f5e5a', ls='', markersize=5,
                          markerfacecolor='none', label='Lower-ranked models'))
    fig.tight_layout()
    bottom = figure_style.stack_legends(fig, [('', handles, 2)], left=0.02)
    fig.subplots_adjust(bottom=(bottom + 0.45) / fig.get_figheight())
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(f'{out}.svg')
    fig.savefig(f'{out}.png', dpi=300)


def build(jobs_file, results, out, table):
    models = scores(pd.read_csv(jobs_file), results)
    draw(models, out)
    models.to_csv(table, index=False)
    best = models[models['model'] == 0]
    print(best[['pair', 'kind', 'partner', 'iptm', 'ipsae', 'ptm']].to_string(index=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--jobs', default='paper/alphafold/jobs.csv')
    parser.add_argument('--results', default=alphafold_results.RESULTS,
                        help='where the server downloads are unpacked, as folds_*')
    parser.add_argument('--table', default='paper/tables/alphafold.csv')
    parser.add_argument('--figure', default='paper/figures/alphafold',
                        help='output path without extension')
    args = parser.parse_args()
    build(args.jobs, args.results, args.figure, args.table)
