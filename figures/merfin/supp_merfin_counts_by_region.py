#!/usr/bin/env python3
"""supp_merfin_counts_by_region.py - supplementary figure: GATK-only and PPR-only variants by region, all and with k-mer support (median
per individual). Output: figures/SuppFig_merfin_counts_by_region.{png,pdf}.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure5_panels as f5  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

f5.style()
R, NV, HG, HR = f5.load()
fig = plt.figure(figsize=(7.1, 3.1))
handles = f5.draw_counts_by_region(fig, NV, R)
fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=6.8)
f5.supp_title(fig, 'Variants found by one pipeline and not by the other, by region, all and with k-mer support')
f5.save(fig, 'SuppFig_merfin_counts_by_region')
print('SUPP_MF_COUNTS_BY_REGION_DONE')
