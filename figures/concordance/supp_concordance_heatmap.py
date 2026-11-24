#!/usr/bin/env python3
"""Supplementary Figure S8: mean genotype concordance of the nine call-set pairs in the genome and the 10 regions, as a heatmap with
the value in each cell, for all variants, SNVs and indels.
Usage: python3 supp_concordance_heatmap.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure3_panels as f3  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

PAIRS = [p for p, _ in f3.PAIR_ROWS if p is not None]
GROUPS = [('Illumina', PAIRS[0:3]), ('Aviti', PAIRS[3:6]), ('Illumina vs Aviti', PAIRS[6:9])]
LABEL = dict((p, lab) for p, lab in f3.PAIR_ROWS if p is not None)
VMIN, VMAX = 30.0, 100.0
f3.style()
D, T = f3.load()
cmap = plt.get_cmap('viridis')
fig, axes = plt.subplots(3, 1, figsize=(7.1, 7.6))
for ax, (typ, tl) in zip(axes, f3.TYPES):
    M = np.array([[100 * D[('genotype', p, r, typ)].mean() for r in f3.ROWS] for p in PAIRS])
    if M.min() < VMIN or M.max() > VMAX:
        raise SystemExit(f'{typ}: values outside {VMIN}-{VMAX}')
    im = ax.imshow(M, cmap=cmap, vmin=VMIN, vmax=VMAX, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            rgb = cmap((M[i, j] - VMIN) / (VMAX - VMIN))[:3]
            light = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] > 0.55
            ax.text(j, i, f'{M[i, j]:.1f}', ha='center', va='center', fontsize=5.8, color='#1f2933' if light else 'white')
    ax.set_xticks(range(len(f3.ROWS)))
    ax.set_xticklabels([f3.ROW_LABEL[r] for r in f3.ROWS])
    ax.xaxis.tick_top()
    ax.tick_params(axis='both', length=0)
    ax.set_yticks(range(len(PAIRS)))
    ax.set_yticklabels([f'{g}: {LABEL[p]}' for g, ps in GROUPS for p in ps], fontsize=6.3)
    for y in (2.5, 5.5):
        ax.axhline(y, color='white', lw=2.5)
    ax.axvline(0.5, color='white', lw=2.5)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(tl, fontweight='bold', loc='left', pad=14)
fig.subplots_adjust(hspace=0.35, right=0.88)
cax = fig.add_axes([0.905, 0.25, 0.018, 0.5])
cb = fig.colorbar(im, cax=cax)
cb.set_label('Genotype concordance (%)', fontsize=7)
cb.ax.tick_params(labelsize=6.5)
f3.save(fig, 'SuppFig_concordance_heatmap')
print('SUPP_CONCORDANCE_HEATMAP_DONE')
