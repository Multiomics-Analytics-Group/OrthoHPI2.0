'''
The confidence score of predicted links against the AlphaFold 3 Server's confidence in
their interface (scripts/build_alphafold_benchmark_jobs.py): the ipTM (A) and ipSAE (B)
of the best-ranked model of each link, coloured by the kind of binding BioGRID shows for
its orthologous groups in human (scripts/build_benchmark_evidence.py) or, with
--colour-by evidence, by the evidence its score mostly rests on, with the random host
proteins' range where their models are in. Runs on whatever results are in so far.
Writes paper/figures/alphafold_benchmark.svg/.png (or alphafold_benchmark_by_evidence)
and paper/tables/alphafold_benchmark.csv.
'''
import argparse
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts import alphafold_results, figure_style
from scripts.build_alphafold_figure import CONFIDENT, FAILED

EVIDENCE = ['experimental', 'databases']
EVIDENCE_LABELS = {'experimental': 'Mostly experimental evidence',
                   'databases': 'Mostly database evidence'}
EVIDENCE_COLORS = {'experimental': '#2a78d6', 'databases': '#eb6834'}
EVIDENCE_MARKERS = {'experimental': 'o', 'databases': 's'}
# drawn in this order, the links without human evidence first, behind the rest
BINDING = ['none found', 'proximity', 'co-complex', 'direct']
BINDING_LABELS = {'direct': 'Direct binding in human', 'co-complex': 'Co-complex only',
                  'proximity': 'Proximity labelling only', 'none found': 'No human experiment'}
BINDING_COLORS = {'direct': '#2a78d6', 'co-complex': '#eb6834', 'proximity': '#4a3aa7',
                  'none found': '#c3c2b7'}
CONTROL_COLOR = '#1baf7a'
METRICS = [('iptm', 'Interface confidence (ipTM)'), ('ipsae', 'ipSAE')]


def scores(jobs, results):
    '''The best-ranked model of every job run so far, with its link.'''
    rows = []
    for job in jobs.to_dict('records'):
        folder = alphafold_results.job_folder(job['job'], results)
        if folder is not None:
            rows.append({**job, **alphafold_results.models(folder)[0]})

    return pd.DataFrame(rows)


def draw(best, column, out):
    figure_style.apply()
    if column == 'evidence':
        groups, colors, labels = EVIDENCE, EVIDENCE_COLORS, EVIDENCE_LABELS
    else:
        groups, colors, labels = BINDING, BINDING_COLORS, BINDING_LABELS
    predicted = best[best['kind'] == 'predicted']
    control = best[best['kind'] == 'control']
    fig, axes = plt.subplots(1, 2, figsize=(figure_style.WIDTH, 3.2), sharex=True)
    for ax, (metric, label), letter in zip(axes, METRICS, 'AB'):
        if metric == 'iptm':
            ax.axhspan(0, FAILED, color='#f0efeb', zorder=0, lw=0)
            ax.axhline(CONFIDENT, color='#b4b2a9', lw=0.6, ls='--', zorder=0)
        if len(control):
            low, high = control[metric].quantile([0.05, 0.95])
            ax.axhspan(low, high, color=CONTROL_COLOR, alpha=0.15, lw=0, zorder=1)
            ax.axhline(control[metric].median(), color=CONTROL_COLOR, lw=0.8, zorder=1)
        notes = []
        for group in groups:
            for evidence in EVIDENCE:
                rows = predicted[(predicted[column] == group) & (predicted['evidence'] == evidence)]
                ax.scatter(rows['score'], rows[metric], s=22, marker=EVIDENCE_MARKERS[evidence],
                           color=colors[group], edgecolors='white', lw=0.5, zorder=3)
            rows = predicted[predicted[column] == group]
            if column == 'evidence' and len(rows) > 2:
                rho, p = spearmanr(rows['score'], rows[metric])
                notes.append(f'{group}: ρ = {rho:.2f} (p = {p:.2g}, n = {len(rows)})')
        if column != 'evidence':
            rho, p = spearmanr(predicted['score'], predicted[metric])
            notes.append(f'ρ = {rho:.2f} (p = {p:.2g}, n = {len(predicted)})')
        ax.text(0.02, 0.98, '\n'.join(notes), transform=ax.transAxes, va='top',
                fontsize=figure_style.FONT, color='#5f5e5a')
        ax.set_ylim(0, 1)
        ax.set_xlim(0.33, 1)
        ax.set_xlabel('Confidence score of the predicted link')
        ax.set_ylabel(label)
        ax.set_title(letter, loc='left', fontweight='bold')
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)

    predicted_counts = predicted[column].value_counts()
    handles = [Line2D([], [], marker='o', color=colors[g], ls='', markersize=5,
                      label=f'{labels[g]} (n = {predicted_counts.get(g, 0)})')
               for g in reversed(groups)]
    if column != 'evidence':
        handles += [Line2D([], [], marker=EVIDENCE_MARKERS[e], color='#5f5e5a', ls='',
                           markersize=5, label=EVIDENCE_LABELS[e]) for e in EVIDENCE]
    if len(control):
        handles.append(Patch(color=CONTROL_COLOR, alpha=0.3,
                             label=f'Random host proteins, 5–95% and median (n = {len(control)})'))
    fig.tight_layout()
    bottom = figure_style.stack_legends(fig, [('', handles, 2)], left=0.06)
    fig.subplots_adjust(bottom=(bottom + 0.5) / fig.get_figheight())
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(f'{out}.svg')
    fig.savefig(f'{out}.png', dpi=300)


def build(jobs_file, results, colour_by, binding_file, out, table):
    jobs = pd.read_csv(jobs_file)
    best = scores(jobs, results)
    if best.empty:
        raise SystemExit(f'no benchmark results in {results}')
    if colour_by == 'binding':
        if not os.path.exists(binding_file):
            raise SystemExit(f'{binding_file} not found: run scripts/build_benchmark_evidence.py')
        binding = pd.read_csv(binding_file)[['link', 'kind']].rename(columns={'kind': 'binding'})
        best = best.merge(binding, on='link', how='left')
    else:
        out = f'{out}_by_evidence'
    draw(best, colour_by, out)
    best.to_csv(table, index=False)
    print(f"{len(best)} of {len(jobs)} jobs run: "
          f"{(best['kind'] == 'predicted').sum()} links, {(best['kind'] == 'control').sum()} controls")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--jobs', default='paper/alphafold/benchmark.csv')
    parser.add_argument('--results', default=alphafold_results.RESULTS,
                        help='where the server downloads are unpacked, as folds_*')
    parser.add_argument('--colour-by', choices=['binding', 'evidence'], default='binding')
    parser.add_argument('--binding', default='paper/tables/alphafold_benchmark_evidence.csv',
                        help='the BioGRID labels of scripts/build_benchmark_evidence.py')
    parser.add_argument('--table', default='paper/tables/alphafold_benchmark.csv')
    parser.add_argument('--figure', default='paper/figures/alphafold_benchmark',
                        help='output path without extension')
    args = parser.parse_args()
    build(args.jobs, args.results, args.colour_by, args.binding, args.figure, args.table)
