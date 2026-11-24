#!/usr/bin/env python3
"""supp_merfin_hg002_giab_classes.py - supplementary figure: k-mer support of the GIAB true positives and of the false positives with a
wrong allele or a wrong genotype on HG002, per call set and platform. Output: figures/SuppFig_merfin_hg002_giab_classes.{png,pdf}.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure5_panels as f5  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

f5.style()
R, NV, HG, HR = f5.load()
fig, axes = plt.subplots(2, 1, figsize=(7.1, 4.6), sharex=True)
handles = f5.draw_support_by_class(axes, HG)
fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.01), ncol=3, fontsize=6.5)
fig.subplots_adjust(hspace=0.12)
f5.supp_title(fig, 'K-mer support of the GIAB true and false positives on HG002')
f5.save(fig, 'SuppFig_merfin_hg002_giab_classes')
print('SUPP_MF_HG002_CLASSES_DONE')
