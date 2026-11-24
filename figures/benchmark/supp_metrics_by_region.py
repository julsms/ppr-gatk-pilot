#!/usr/bin/env python3
"""Supplementary Figures S1 (PASS) and S2 (--filter ALL: GATK before VQSR; the DeepVariant pipelines have no filter step, so their PASS
values are used): precision, recall and F1 of GATK, BWA-MEM + DeepVariant and PPR on Illumina and Aviti in the 10 regions, as dots.
Inputs: benchmark_data.tsv. Outputs: figures/SuppFig_metrics_by_region_<PASS|ALL>.png (600 dpi) and .pdf.
Usage: python3 supp_metrics_by_region.py [--filter PASS|ALL]
"""
import csv
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'benchmark_data.tsv')
REGIONS = ['BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']
REGION_LABEL = dict({r: r for r in REGIONS}, BB='BB (all)')
PIPES = ['GATK', 'BWA-MEM + DeepVariant', 'PPR']
COLOR = {'GATK': '#882255', 'BWA-MEM + DeepVariant': '#332288', 'PPR': '#117733'}   # Tol muted, as Figure 1
MARKER = {'Illumina': 'o', 'Aviti': 's'}
TYPES = [('SNP', 'SNVs'), ('INDEL', 'Indels')]
METRICS = [('precision', 'Precision'), ('recall', 'Recall'), ('f1', 'F1')]
INK, GRID = '#1f2933', '#e3e7eb'

filt = 'PASS'
if '--filter' in sys.argv:
    filt = sys.argv[sys.argv.index('--filter') + 1]
if filt not in ('PASS', 'ALL'):
    sys.exit('--filter must be PASS or ALL')

rows = list(csv.DictReader(open(DATA), delimiter='\t'))
if len(rows) != 160 or sorted({r['region'] for r in rows}) != sorted(REGIONS):
    sys.exit(f'unexpected table {DATA}: {len(rows)} rows')
val = {}
for r in rows:
    val[(r['platform'], r['pipeline'], r['region'], r['type'], r['filter'])] = {m: float(r[m]) for m, _ in METRICS}


def value(plat, pipe, reg, typ, metric):
    f = filt if pipe == 'GATK' else 'PASS'   # the DeepVariant pipelines have no filter step
    return val[(plat, pipe, reg, typ, f)][metric]


plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10, 'axes.edgecolor': '#9aa5b1',
                     'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'axes.spines.top': False,
                     'axes.spines.right': False, 'legend.frameon': False, 'axes.unicode_minus': False})
SERIES = [(plat, pipe) for plat in ('Illumina', 'Aviti') for pipe in PIPES]
fig, axes = plt.subplots(2, 3, figsize=(11, 8.2), sharey=True)
offs = np.linspace(-0.32, 0.32, len(SERIES))
ys = np.arange(len(REGIONS))[::-1]
for i, (typ, typ_label) in enumerate(TYPES):
    for j, (metric, metric_label) in enumerate(METRICS):
        ax = axes[i, j]
        vals = []
        for k, (plat, pipe) in enumerate(SERIES):
            x = [value(plat, pipe, reg, typ, metric) for reg in REGIONS]
            vals += x
            ax.plot(x, ys + offs[k], ls='', marker=MARKER[plat], ms=4.2, mew=1.1, color=COLOR[pipe],
                    markerfacecolor=COLOR[pipe] if plat == 'Illumina' else 'white')
        for y in ys:
            ax.axhspan(y - 0.45, y + 0.45, color=GRID if y % 2 == 0 else 'white', zorder=0, lw=0)
        lo = min(vals)
        ax.set_xlim(max(0.0, lo - 0.03 * (1 - lo + 0.05)), 1.0005)
        ax.set_ylim(-0.6, len(REGIONS) - 0.4)
        ax.set_yticks(ys)
        if j == 0:
            ax.set_yticklabels([REGION_LABEL[r] for r in REGIONS])
            ax.set_ylabel(typ_label, fontweight='bold')
        else:
            ax.tick_params(labelleft=False)
        ax.tick_params(axis='y', length=0)
        ax.grid(axis='x', color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        if i == 0:
            ax.set_title(metric_label, fontweight='bold')
gatk_label = 'GATK' if filt == 'PASS' else 'GATK before VQSR'
handles = [Line2D([], [], ls='', marker=MARKER[plat], ms=6, mew=1.3, color=COLOR[pipe],
                  markerfacecolor=COLOR[pipe] if plat == 'Illumina' else 'white',
                  label=f'{gatk_label if pipe == "GATK" else pipe}, {plat}') for pipe in PIPES for plat in ('Illumina', 'Aviti')]
fig.legend(handles=handles, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.05), handletextpad=0.3, columnspacing=1.2, fontsize=8)
fig.subplots_adjust(bottom=0.1, wspace=0.08, hspace=0.12)
os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True)
for ext, kw in (('png', {'dpi': 600}), ('pdf', {})):
    out = os.path.join(HERE, 'figures', f'SuppFig_metrics_by_region_{filt}.{ext}')
    fig.savefig(out, bbox_inches='tight', facecolor='white', **kw)
    print('written', out)
print('SUPP_METRICS_DONE')
