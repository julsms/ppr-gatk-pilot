#!/usr/bin/env python3
"""supp_uniquevariants_singletons.py - Supplementary Figure S17: the share of singletons among the unique variants (carried by no
other individual in the calls of that pipeline; 22 individuals, one of the two first-degree relatives left out), in the genome and
the 10 regions; all variants. Input: uv_singletons.tsv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure4_panels as f4  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

SETS = [('GATK', 'PPR', ''), ('PPR', 'GATK', ''), ('DV', 'PPR', '////'), ('PPR', 'DV', '////'), ('PPR', 'BOTH', 'xxxx')]
f4.style()
S = {(r['set'], r['platform'], r['region'], r['type']): r['share'] for r in f4.read('uv_singletons.tsv')}
fig, axes = plt.subplots(2, 1, figsize=(7.1, 4.2), sharex=True, sharey=True)
wb = 0.8 / len(SETS)
for pi, (plat, pl) in enumerate(f4.PLATFORMS):
    ax = axes[pi]
    for g, r in enumerate(f4.ROWS):
        for s, (a, b, hatch) in enumerate(SETS):
            v = S[(f'{a}vs{b}', plat, r, 'ALL')]
            if v == 'NA':
                continue
            ax.bar(g - 0.4 + wb * (s + 0.5), 100 * float(v), width=wb * 0.9, color=f4.SHADE[a][pi], hatch=hatch, lw=0.4,
                   edgecolor=(f4.SHADE[a][0] if pi else 'white') if hatch else 'white')
    ax.set_ylabel(f'Singletons (%)\n{pl}')
    ax.grid(axis='y', color=f4.GRID, lw=0.6)
    ax.set_axisbelow(True)
axes[-1].set_xticks(range(len(f4.ROWS)))
axes[-1].set_xticklabels([f4.LAB[r] for r in f4.ROWS])
# two legend columns, filled column by column: Illumina on the left, Aviti on the right, one set per row
handles = [Patch(facecolor=f4.SHADE[a][pi], edgecolor=(f4.SHADE[a][0] if pi else 'white') if hatch else 'white', hatch=hatch,
                 label=f'{f4.set_label(a, b)}, {pl}') for pi, (_, pl) in enumerate(f4.PLATFORMS) for a, b, hatch in SETS]
fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.03), ncol=2, fontsize=6.2)
fig.subplots_adjust(hspace=0.12)
f4.supp_title(fig, 'Singletons among the unique variants, by region')
f4.save(fig, 'SuppFig_uniquevariants_singletons')
print('SUPP_UV_SINGLETONS_DONE')
