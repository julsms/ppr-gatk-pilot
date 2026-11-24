#!/usr/bin/env python3
"""supp_uniquevariants_hg002_precision.py - Supplementary Figure S4: the precision of the unique calls of HG002 against GIAB v5.0q,
inside the benchmark regions (BB), all variants, SNVs and indels. Precision = TP / (TP + FP + FP_gt): a call with the right allele
and a wrong genotype counts as false. Input: uv_hg002_precision.tsv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure4_panels as f4  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

SETS = [('GATK', 'PPR', ''), ('PPR', 'GATK', ''), ('DV', 'PPR', '////'), ('PPR', 'DV', '////'), ('PPR', 'BOTH', 'xxxx')]
f4.style()
P = {(r['set'], r['platform'], r['region'], r['type']): r['precision'] for r in f4.read('uv_hg002_precision.tsv')}
fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.8), sharey=True)
for ax, (typ, tl) in zip(axes, f4.TYPES):
    for g, (a, b, hatch) in enumerate(SETS):
        for pi, (plat, _) in enumerate(f4.PLATFORMS):
            v = 100 * float(P[(f'{a}vs{b}', plat, 'genome', typ)])
            x = g + (pi - 0.5) * 0.38
            ax.bar(x, v, width=0.36, color=f4.SHADE[a][pi], hatch=hatch, lw=0.4,
                   edgecolor=(f4.SHADE[a][0] if pi else 'white') if hatch else 'white')
            ax.text(x, v + 1, f'{v:.0f}', ha='center', va='bottom', fontsize=5.6, color=f4.INK)
    ax.set_xticks([])
    ax.set_ylim(0, 105)
    ax.set_title(tl, fontweight='bold')
    ax.grid(axis='y', color=f4.GRID, lw=0.6)
    ax.set_axisbelow(True)
axes[0].set_ylabel('Precision on HG002, inside BB (%)')
# two legend columns, filled column by column: Illumina on the left, Aviti on the right, one set per row
handles = [Patch(facecolor=f4.SHADE[a][pi], edgecolor=(f4.SHADE[a][0] if pi else 'white') if hatch else 'white', hatch=hatch,
                 label=f'{f4.set_label(a, b)}, {pl}') for pi, (_, pl) in enumerate(f4.PLATFORMS) for a, b, hatch in SETS]
fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.03), ncol=2, fontsize=6.0)
fig.subplots_adjust(wspace=0.06)
f4.supp_title(fig, 'Precision of the unique calls on HG002, inside the benchmark regions')
f4.save(fig, 'SuppFig_uniquevariants_hg002_precision')
print('SUPP_UV_HG002_DONE')
