#!/usr/bin/env python3
"""supp_genomescope_het_vs_length.py - supplementary figure: heterozygosity against the genome length estimate, one point per
individual and platform, Spearman rho per platform. Output: figures/SuppFig_genomescope_het_vs_length.{pdf,png}.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genomescope_panels as gp  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

gp.style()
V, N, ST, COR = gp.load()
fig, ax = plt.subplots(figsize=(3.4, 2.6))
gp.draw_het_vs_length(ax, V, N, COR)
gp.save(fig, 'SuppFig_genomescope_het_vs_length')
print('SUPP_GENOMESCOPE_HET_VS_LENGTH_DONE')
