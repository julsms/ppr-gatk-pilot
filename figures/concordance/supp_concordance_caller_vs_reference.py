#!/usr/bin/env python3
"""Supplementary Figure S10: on each platform, genotype concordance of GATK vs DV (the caller changes) and of DV vs PPR (the reference
changes) by region, with the paired test of the two.
Usage: python3 supp_concordance_caller_vs_reference.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure3_panels as f3  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

PLATFORMS = [('illumina', 'Illumina'), ('aviti', 'Aviti')]
f3.style()
D, T = f3.load()
ypos = f3.region_rows()
fig, axes = plt.subplots(2, 3, figsize=(7.1, 5.9), sharey=True)
for i, (plat, pl) in enumerate(PLATFORMS):
    for j, (typ, tl) in enumerate(f3.TYPES):
        ax = axes[i, j]
        for r in f3.ROWS:
            a = 100 * D[('genotype', f'GATK_DV_{plat}', r, typ)].mean()
            b = 100 * D[('genotype', f'DV_PPR_{plat}', r, typ)].mean()
            ax.plot([a, b], [ypos[r]] * 2, color='#9aa5b1', lw=1.1, zorder=1)
            ax.plot(a, ypos[r], 'o', ms=4.6, color=f3.STEP_COLOR['GATK_DV'], zorder=3)
            ax.plot(b, ypos[r], 'o', ms=4.6, color=f3.STEP_COLOR['DV_PPR'], zorder=3)
        f3.row_axes(ax, ypos, j == 0)
        f3.star_column(ax, T, f'H3_PPRDV_vs_DVGATK_{plat}', typ, ypos)
        ax.set_title(f'{pl}, {tl}', fontweight='bold')
fig.supxlabel('Genotype concordance on one platform (%), mean of the 23 individuals', fontsize=7.5, y=0.05)
h = [Line2D([], [], ls='', marker='o', ms=5, color=f3.STEP_COLOR[s], label=f3.STEP_LABEL[s]) for s in ('GATK_DV', 'DV_PPR')]
fig.legend(handles=h, loc='upper center', bbox_to_anchor=(0.5, 0.035), ncol=2, fontsize=6.4)
fig.subplots_adjust(wspace=0.38, hspace=0.3, bottom=0.11)
os.makedirs(os.path.join(f3.HERE, 'figures'), exist_ok=True)
for ext, kw in (('png', {'dpi': 600}), ('pdf', {})):
    fig.savefig(os.path.join(f3.HERE, 'figures', f'SuppFig_concordance_caller_vs_reference.{ext}'), bbox_inches='tight',
                facecolor='white', **kw)
print('written figures/SuppFig_concordance_caller_vs_reference')
print('SUPP_CONCORDANCE_S9_DONE')
