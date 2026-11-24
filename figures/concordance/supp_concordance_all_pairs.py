#!/usr/bin/env python3
"""Supplementary Figure S9: genotype concordance of every call-set pair in every region, mean of the 23 individuals.
Usage: python3 supp_concordance_all_pairs.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure3_panels as f3  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ORDER = ([(f'{s}_{p}', s, p) for s in ('GATK_DV', 'DV_PPR', 'GATK_PPR') for p in ('illumina', 'aviti')]
         + [(f'{s}_illumina_aviti', s, 'cross') for s in f3.PIPES])
f3.style()
D, T = f3.load()
ypos = f3.region_rows()
offs = np.linspace(0.33, -0.33, len(ORDER))
fig, axes = plt.subplots(1, 3, figsize=(7.1, 5.0), sharey=True)
for ax, (typ, tl) in zip(axes, f3.TYPES):
    for (pair, s, p), o in zip(ORDER, offs):
        for r in f3.ROWS:
            v = 100 * D[('genotype', pair, r, typ)].mean()
            if p == 'cross':
                ax.plot(v, ypos[r] + o, 'D', ms=3.2, color=f3.PIPE_COLOR[s])
            else:
                ax.plot(v, ypos[r] + o, 'o' if p == 'illumina' else 's', ms=3.2, color=f3.STEP_COLOR[s],
                        mfc=f3.STEP_COLOR[s] if p == 'illumina' else 'white', mew=0.8)
    f3.row_axes(ax, ypos, ax is axes[0])
    ax.set_title(tl, fontweight='bold')
for ax in axes:
    ax.xaxis.set_major_locator(f3.MaxNLocator(nbins=4, integer=True))
fig.supxlabel('Genotype concordance (%), mean of the 23 individuals', fontsize=7.5, y=0.03)
fig.legend(handles=f3.legend_pairs()[:-1], loc='upper center', bbox_to_anchor=(0.5, 0.005), ncol=4, fontsize=6.4)
fig.subplots_adjust(wspace=0.14, bottom=0.12)
os.makedirs(os.path.join(f3.HERE, 'figures'), exist_ok=True)
for ext, kw in (('png', {'dpi': 600}), ('pdf', {})):
    fig.savefig(os.path.join(f3.HERE, 'figures', f'SuppFig_concordance_all_pairs.{ext}'), bbox_inches='tight', facecolor='white', **kw)
print('written figures/SuppFig_concordance_all_pairs')
print('SUPP_CONCORDANCE_S8_DONE')
