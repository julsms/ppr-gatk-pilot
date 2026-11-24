#!/usr/bin/env python3
"""supp_uniquevariants_1000g.py - Supplementary Figure S14: the allele frequency in 1000 Genomes (phase 3, global) of the
pipeline-unique variants, as Figure 4b for gnomAD genomes; the share found in 1000 Genomes is written in each panel.
Input: uv_af_hist.tsv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure4_panels as f4  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

f4.style()
C, H, Q = f4.load()
fig, axes = plt.subplots(2, 2, figsize=(7.1, 4.4), sharex=True)
for i, pair in enumerate(f4.PAIRS):
    for j, (plat, pl) in enumerate(f4.PLATFORMS):
        f4.draw_af(axes[i][j], H, pair, plat, source='1000G', ylabel=(j == 0))
        axes[i][j].set_title(f'{pair[0]} vs {pair[1]}, {pl}', fontweight='bold')
fig.supxlabel('Allele frequency in 1000 Genomes', fontsize=7, y=0.02)
fig.subplots_adjust(hspace=0.32, wspace=0.2)
f4.supp_title(fig, 'Allele frequency in 1000 Genomes of the pipeline-unique variants')
f4.save(fig, 'SuppFig_uniquevariants_1000g')
print('SUPP_UV_1000G_DONE')
