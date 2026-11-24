#!/usr/bin/env python3
"""supp_merfin_regions_cohort.py - supplementary figure: k-mer support by region in the cohort, mean of the 23 individuals, all variants,
Illumina and Aviti. Output: figures/SuppFig_merfin_regions_cohort.{png,pdf}.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure5_panels as f5  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

f5.style()
R, NV, HG, HR = f5.load()
fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.7), sharey=True)
im = f5.draw_region_heatmaps(axes, R)
cb = fig.colorbar(im, ax=axes, fraction=0.03, pad=0.02)
cb.set_label('Records with k-mer support (%), mean of the individuals', fontsize=6.3)
cb.ax.tick_params(labelsize=6)
f5.supp_title(fig, 'K-mer support by region in the cohort')
f5.save(fig, 'SuppFig_merfin_regions_cohort')
print('SUPP_MF_REGIONS_COHORT_DONE')
