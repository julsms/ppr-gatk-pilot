#!/usr/bin/env python3
"""supp_merfin_per_individual.py - supplementary figure: k-mer support of each individual, genome, all variants (violins: left half
Illumina, right half Aviti). Output: figures/SuppFig_merfin_per_individual.{png,pdf}.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure5_panels as f5  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

f5.style()
R, NV, HG, HR = f5.load()
fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.9), gridspec_kw={'width_ratios': [3, 5.2], 'wspace': 0.25})
handles = f5.draw_violins(axes, R)
fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=6, fontsize=6.3)
f5.supp_title(fig, 'K-mer support per individual')
f5.save(fig, 'SuppFig_merfin_per_individual')
print('SUPP_MF_PER_INDIVIDUAL_DONE')
