#!/usr/bin/env python3
"""supp_uniquevariants_why_gatk_missed.py - Supplementary Figure S15: the GATK status of the PPR-only variants (vs GATK), by region:
removed by VQSR, by the ExcessHet filter or by another filter; called by GATK in another individual and 0/0 (or missing) in this
one; not called by GATK in any individual. Colours of Paul Tol's muted scheme. Input: uv_r119.tsv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure4_panels as f4  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

CLASSES = [('vqsr', 'GATK call removed by VQSR', '#882255'), ('excesshet', 'GATK call removed by the ExcessHet filter', '#CC6677'),
           ('other', 'other filter or not matched', '#DDDDDD'),
           ('ref_in_person', 'allele called by GATK in another individual, this one 0/0', '#88CCEE'),
           ('missing_in_person', 'allele called by GATK in another individual, this one missing', '#332288'),
           ('not_called', 'allele not called by GATK in any individual', '#999933')]
f4.style()
share = {(r['platform'], r['region'], r['type'], r['class']): float(r['share']) for r in f4.read('uv_r119.tsv')}
fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.4), sharey=True)
yy = {r: -i - (0.6 if i else 0.0) for i, r in enumerate(f4.ROWS)}
for ax, (plat, pl) in zip(axes, f4.PLATFORMS):
    for r in f4.ROWS:
        left = 0.0
        for key, _, c in CLASSES:
            w = 100 * share[(plat, r, 'ALL', key)]
            ax.barh(yy[r], w, left=left, height=0.75, color=c, edgecolor='white', lw=0.4)
            left += w
        if abs(left - 100) > 1e-3:
            raise SystemExit(f'{plat} {r}: the classes add up to {left:.4f} %')
    ax.set_yticks([yy[r] for r in f4.ROWS])
    ax.set_yticklabels([f4.LAB[r] for r in f4.ROWS])
    ax.tick_params(axis='y', length=0)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_title(pl, fontweight='bold')
    ax.set_xlabel('Share of the PPR-only variants (vs GATK) (%)')
fig.legend(handles=[Patch(color=c, label=lab) for _, lab, c in CLASSES], loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=2,
           fontsize=6.3)
fig.subplots_adjust(wspace=0.12)
f4.supp_title(fig, 'GATK status of the PPR-only variants, by region')
f4.save(fig, 'SuppFig_uniquevariants_why_gatk_missed')
print('SUPP_UV_WHY_GATK_DONE')
