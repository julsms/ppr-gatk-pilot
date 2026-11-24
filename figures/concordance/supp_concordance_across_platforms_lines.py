#!/usr/bin/env python3
"""Supplementary Figure S11: genotype concordance across platforms of each pipeline, all variants, by region (mean and 25th to 75th
percentile), with the tests of PPR vs GATK and of PPR vs DV.
Usage: python3 supp_concordance_across_platforms_lines.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure3_panels as f3  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.transforms as mtransforms  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

TYPE = 'ALL'
DODGE = {'GATK': -0.14, 'DV': 0.0, 'PPR': 0.14}
f3.style()
D, T = f3.load()
xs = np.array([-1.3] + list(range(len(f3.REGIONS))), dtype=float)   # the genome apart from the regions
fig, ax = plt.subplots(figsize=(7.1, 3.3))
for s in f3.PIPES:
    vals = [100 * D[('genotype', f'{s}_illumina_aviti', r, TYPE)] for r in f3.ROWS]
    mean = np.array([v.mean() for v in vals])
    q1 = np.array([np.percentile(v, 25) for v in vals])
    q3 = np.array([np.percentile(v, 75) for v in vals])
    x = xs + DODGE[s]
    ax.vlines(x, q1, q3, color=f3.PIPE_COLOR[s], lw=2.2, alpha=0.35, zorder=2)
    ax.plot(x[1:], mean[1:], color=f3.PIPE_COLOR[s], lw=1.6, zorder=3)
    ax.plot(x, mean, ls='', marker='D', ms=4.6, color=f3.PIPE_COLOR[s], mec='white', mew=0.6, zorder=4)
ax.axvline(-0.65, color='#c3cbd3', lw=0.8)
trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
for (cid, other), y in zip(f3.H2_MARKS, (1.075, 1.01)):
    for x, r in zip(xs, f3.ROWS):
        p = float(T[(cid, 'genotype', r, TYPE)]['p_holm'])
        ax.text(x, y, f3.stars(p), transform=trans, ha='center', va='bottom', fontsize=7,
                fontweight='bold' if p < 0.05 else 'normal', color=f3.PIPE_COLOR[other])
ax.set_xticks(xs)
ax.set_xticklabels([f3.ROW_LABEL[r] for r in f3.ROWS], fontsize=7)
ax.set_xlim(-1.9, len(f3.REGIONS) - 0.5)
ax.set_ylabel('Genotype concordance across platforms (%)')
ax.yaxis.set_major_locator(f3.MaxNLocator(nbins=6, integer=True))
ax.grid(axis='y', color=f3.GRID, lw=0.6)
ax.set_axisbelow(True)
pipes = [Line2D([], [], lw=1.6, marker='D', ms=4.6, color=f3.PIPE_COLOR[s], mec='white', mew=0.6, label=f'{s}: Illumina vs Aviti')
         for s in f3.PIPES]
iqr = Line2D([], [], lw=2.2, alpha=0.35, color='#555555', label='25th to 75th percentile of the 23 individuals')
# one row: the three pipelines and the percentile bars (the marks are explained in the figure legend of the manuscript)
fig.legend(handles=pipes + [iqr], loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=4, fontsize=6.6)
os.makedirs(os.path.join(f3.HERE, 'figures'), exist_ok=True)
for ext, kw in (('png', {'dpi': 600}), ('pdf', {})):
    fig.savefig(os.path.join(f3.HERE, 'figures', f'SuppFig_concordance_across_platforms_lines.{ext}'), bbox_inches='tight',
                facecolor='white', **kw)
print('written figures/SuppFig_concordance_across_platforms_lines')
print('SUPP_CONCORDANCE_S10_DONE')
