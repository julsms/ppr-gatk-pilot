#!/usr/bin/env python3
"""supp_genomescope_estimates.py - supplementary figure: the four GenomeScope 2.0 estimates on Illumina and Aviti, 21 individuals
(paired lines, medians, Holm-corrected paired Wilcoxon tests). Output: figures/SuppFig_genomescope_estimates.{pdf,png}.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genomescope_panels as gp  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

gp.style()
V, N, ST, COR = gp.load()
fig, axes = plt.subplots(1, 4, figsize=(7.1, 2.5))
gp.draw_paired_estimates(axes, V, N, ST)
fig.subplots_adjust(wspace=0.55)
gp.save(fig, 'SuppFig_genomescope_estimates')
print('SUPP_GENOMESCOPE_ESTIMATES_DONE')
