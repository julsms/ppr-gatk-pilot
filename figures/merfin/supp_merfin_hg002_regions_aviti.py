#!/usr/bin/env python3
"""supp_merfin_hg002_regions_aviti.py - supplementary figure: Figure 5d for Aviti (k-mer support and GIAB precision by region on HG002,
PPR-only and GATK-only variants). Output: figures/SuppFig_merfin_hg002_regions_aviti.{png,pdf}.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure5_panels as f5  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

f5.style()
R, NV, HG, HR = f5.load()
fig, axes = plt.subplots(2, 1, figsize=(7.1, 4.2), sharex=True)
handles = f5.draw_regions_lines(axes, HR, 'aviti')
axes[0].tick_params(axis='x', labelbottom=False)
axes[1].set_ylabel('Records inside the\nGIAB benchmark (%)')
fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=2, fontsize=6.5)
fig.text(0.5, -0.045, 'Grey numbers under the points: share of the records inside the GIAB benchmark.', ha='center', va='top',
         fontsize=6, color=f5.MUTED)
fig.subplots_adjust(hspace=0.35)
f5.supp_title(fig, 'K-mer support and GIAB precision by region on HG002, Aviti')
f5.save(fig, 'SuppFig_merfin_hg002_regions_aviti')
print('SUPP_MF_HG002_REGIONS_AVITI_DONE')
