#!/usr/bin/env python3
"""supp_depth_by_region.py - supplementary figure: Spearman rho of depth with genotype concordance for the genome and the 10 regions
(a, b within Illumina and Aviti; c between platforms, against the lower depth); filled markers for Holm-adjusted p < 0.05.
Output: figures/SuppFig_depth_by_region.{pdf,png}.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure6_panels as f6  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

f6.style()
fig, axes = plt.subplots(3, 1, figsize=(7.087, 6.3), sharex=True)
fig.subplots_adjust(hspace=0.32, left=0.08, right=0.80, top=0.93, bottom=0.06)
for ax, (pl, plname), tag in zip(axes[:2], f6.PLATS, 'ab'):
    rows = [(lab, f6.PAIR_COL[pid], mk, (lambda st, pid=pid, pl=pl: ('within', f'depth {pl}', pid, pl, 'ALL', st)))
            for (pid, lab), mk in zip(f6.PAIRS, f6.MARK)]
    hand = f6.draw_regions(ax, rows, f'{tag}  Concordance within a platform, {plname}')
    ax.legend(handles=hand, loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=6.3)
rows = [(s, f6.PIPE_COL[s], mk, (lambda st, s=s: ('across', 'lower depth', f'{s} Illumina vs Aviti', 'both', 'ALL', st)))
        for s, mk in zip(f6.PIPES, f6.MARK)]
hand = f6.draw_regions(axes[2], rows, 'c  Concordance between Illumina and Aviti, against the lower of the two depths')
axes[2].legend(handles=hand, loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=6.3)
fig.legend(handles=[Line2D([], [], color=f6.INK, marker='o', ls='', ms=4.2, label='Holm-adjusted p < 0.05'),
                    Line2D([], [], color=f6.INK, marker='o', ls='', ms=4.2, mfc='white', label='not significant')],
           loc='upper left', bbox_to_anchor=(0.815, 0.12), fontsize=6.3)
fig.text(0.08, 0.985, 'Correlation of sequencing depth with concordance, by region', ha='left', va='top', fontsize=8.5,
         fontweight='bold', color=f6.INK)
f6.save(fig, 'SuppFig_depth_by_region')
w, h = f6.size_mm(fig)
plt.close(fig)
print(f'written figures/SuppFig_depth_by_region ({w:.0f} x {h:.0f} mm; PNG at {f6.DPI} dpi)')
print('SUPP_DEPTH_BY_REGION_DONE')
