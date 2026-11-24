#!/usr/bin/env python3
"""Supplementary Figure S28: computing performance of GATK and PPR per sample, Illumina and Aviti. a execution time, b CPU time,
c peak memory of the most memory-intensive process, d peak memory summed over processes. Bars: mean; whiskers: SD.
Input: performance_data.tsv (performance_table.py). Output: figures/SuppFig_performance.png (600 dpi) and .pdf.
"""
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SHADE = {'GATK': ('#b5452f', '#e07b91'), 'PPR': ('#2e7d4f', '#7cb342')}   # the colours of Figure 1: (Illumina, Aviti)
PLATFORMS = [('illumina', 'Illumina'), ('aviti', 'Aviti')]
PANELS = [('time_h', 'Execution time (h)', 1), ('cpu_h', 'CPU time (h)', 0),
          ('max_rss_gb', 'Peak memory (GB)', 0), ('total_rss_gb', 'Summed peak memory (GB)', 0)]
INK, GRID = '#1f2933', '#e3e7eb'

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.labelsize': 7.5, 'xtick.labelsize': 7,
                     'ytick.labelsize': 7, 'axes.edgecolor': '#9aa5b1', 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False})
V = defaultdict(list)
with open(os.path.join(HERE, 'performance_data.tsv'), newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        V[(r['pipeline'], r['platform'], r['metric'])].append(float(r['value']))
fig, axes = plt.subplots(2, 2, figsize=(7.1, 4.8))
for ax, (metric, ylabel, dec), tag in zip(axes.flat, PANELS, 'abcd'):
    top = 0.0
    for g, (plat, pl) in enumerate(PLATFORMS):
        for k, pipe in enumerate(('GATK', 'PPR')):
            v = np.array(V[(pipe, plat, metric)])
            m, sd = v.mean(), v.std(ddof=1)
            x = g + (k - 0.5) * 0.36
            ax.bar(x, m, width=0.34, color=SHADE[pipe][0 if plat == 'illumina' else 1], edgecolor='white', lw=0.5)
            ax.errorbar(x, m, yerr=sd, color=INK, lw=0.8, capsize=2.5)
            ax.text(x, m + sd + 0.02 * max(m + sd, 1), f'{m:,.{dec}f}', ha='center', va='bottom', fontsize=6.5, fontweight='bold',
                    color=INK)
            top = max(top, m + sd)
    ax.set_xticks(range(len(PLATFORMS)))
    ax.set_xticklabels([pl for _, pl in PLATFORMS])
    ax.tick_params(axis='x', length=0)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, top * 1.22)
    ax.grid(axis='y', color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.text(-0.2, 1.04, tag, transform=ax.transAxes, fontsize=11, fontweight='bold', va='bottom')
handles = [Patch(color=SHADE[p][i], label=f'{p}, {pl}') for p in ('GATK', 'PPR') for i, (_, pl) in enumerate(PLATFORMS)]
fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.02), ncol=4, fontsize=7)
fig.subplots_adjust(hspace=0.42, wspace=0.3)
os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True)
for ext, kw in (('png', {'dpi': 600}), ('pdf', {})):
    fig.savefig(os.path.join(HERE, 'figures', f'SuppFig_performance.{ext}'), bbox_inches='tight', facecolor='white', **kw)
print('written figures/SuppFig_performance')
print('SUPP_PERFORMANCE_DONE')
