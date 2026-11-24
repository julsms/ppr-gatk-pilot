#!/usr/bin/env python3
"""Supplementary Figures S5 (Illumina) and S6 (Aviti): F1 of PASS calls on HG002 per region for the five set-ups of the graph
decomposition, SNVs left and indels right; names, order and colours from graph_setups.py.
Inputs: graph_arms_data.tsv. Outputs: figures/SuppFig_graph_arms_F1_<Illumina|Aviti>.png (600 dpi) and .pdf.
Usage: python3 supp_graph_arms_f1.py [--platform Illumina|Aviti]
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
sys.path.insert(0, HERE)
import graph_setups as gs  # noqa: E402

DATA = os.path.join(HERE, 'graph_arms_data.tsv')
REGIONS = gs.REGIONS
REGION_LABEL = dict({r: r for r in REGIONS}, BB='BB (all)')
SETUPS, COLOR, LABEL = gs.SETUPS, gs.COLOR, gs.LABEL
TYPES = [('SNP', 'SNVs'), ('INDEL', 'Indels')]
INK, GRID = '#1f2933', '#e3e7eb'

platforms = ['Illumina', 'Aviti']
if '--platform' in sys.argv:
    platforms = [sys.argv[sys.argv.index('--platform') + 1]]
    if platforms[0] not in ('Illumina', 'Aviti'):
        sys.exit('--platform must be Illumina or Aviti')

rows = list(csv.DictReader(open(DATA), delimiter='\t'))
if len(rows) != 200 or sorted({r['region'] for r in rows}) != sorted(REGIONS) or sorted({r['setup'] for r in rows}) != sorted(SETUPS):
    sys.exit(f'unexpected table {DATA}: {len(rows)} rows')
f1 = {(r['platform'], r['setup'], r['region'], r['type']): float(r['f1']) for r in rows}

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10, 'axes.edgecolor': '#9aa5b1',
                     'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'axes.spines.top': False,
                     'axes.spines.right': False, 'legend.frameon': False, 'axes.unicode_minus': False})
os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True)
for plat in platforms:
    fig, axes = plt.subplots(1, 2, figsize=(11, 6.8), sharey=True)
    ys = np.arange(len(REGIONS))[::-1]
    offs = np.linspace(-0.28, 0.28, len(SETUPS))
    for ax, (typ, typ_label) in zip(axes, TYPES):
        vals = []
        for k, s in enumerate(SETUPS):
            x = [f1[(plat, s, r, typ)] for r in REGIONS]
            vals += x
            ax.plot(x, ys + offs[k], ls='', marker='o', ms=5, color=COLOR[s])
        for y in ys:
            ax.axhspan(y - 0.45, y + 0.45, color=GRID if y % 2 == 0 else 'white', zorder=0, lw=0)
        ax.set_xlim(min(vals) - 0.01, 1.0005)
        ax.set_ylim(-0.6, len(REGIONS) - 0.4)
        ax.set_yticks(ys)
        ax.set_yticklabels([REGION_LABEL[r] for r in REGIONS])
        ax.tick_params(axis='y', length=0)
        ax.set_xlabel('F1 (PASS)')
        ax.set_title(typ_label, fontweight='bold')
        ax.grid(axis='x', color=GRID, lw=0.6)
        ax.set_axisbelow(True)
    handles = [Line2D([], [], ls='', marker='o', ms=7, color=COLOR[s], label=LABEL[s]) for s in SETUPS]
    fig.legend(handles=handles, loc='lower left', ncol=1, bbox_to_anchor=(0.1, -0.17), fontsize=8.5,
               title=gs.LEGEND_TITLE, title_fontsize=8.5, alignment='left')
    fig.subplots_adjust(bottom=0.12, wspace=0.06)
    for ext, kw in (('png', {'dpi': 600}), ('pdf', {})):
        out = os.path.join(HERE, 'figures', f'SuppFig_graph_arms_F1_{plat}.{ext}')
        fig.savefig(out, bbox_inches='tight', facecolor='white', **kw)
        print('written', out)
    plt.close(fig)
print('SUPP_GRAPH_ARMS_DONE')
