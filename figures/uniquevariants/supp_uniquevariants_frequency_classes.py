#!/usr/bin/env python3
"""supp_uniquevariants_frequency_classes.py - Supplementary Figure S13: the allele frequency in gnomAD genomes of the pipeline-unique
variants, in the genome and the 10 regions: the share of the distinct alleles of the cohort not in gnomAD genomes, below 0.1 %, 0.1 to
1 %, 1 to 5 %, 5 % or more. Input: uv_af_classes.tsv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure4_panels as f4  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

CLASSES = [('absent', 'not in gnomAD genomes', '#b9bec5'), ('lt0.1pct', 'below 0.1 %', '#2d1e4a'), ('0.1-1pct', '0.1 to 1 %', '#5b3f8c'),
           ('1-5pct', '1 to 5 %', '#9a7cc4'), ('ge5pct', '5 % or more', '#d9cdee')]
SETS = [('GATK', 'PPR'), ('PPR', 'GATK'), ('DV', 'PPR'), ('PPR', 'DV'), ('PPR', 'BOTH')]
f4.style()
n = {}
for r in f4.read('uv_af_classes.tsv'):
    if r['source'] == 'gnomADg' and r['type'] == 'ALL':
        n[(r['set'], r['platform'], r['region'], r['class'])] = int(r['alleles'])
labels, ypos, y = [], [], 0.0
for a, b in SETS:
    for plat, pl in f4.PLATFORMS:
        labels.append(f'{f4.set_label(a, b)}, {pl}')
        ypos.append(y)
        y -= 1.0
    y -= 0.5
fig, axes = plt.subplots(3, 4, figsize=(7.1, 7.4), sharey=True)
for k, ax in enumerate(axes.flat):
    if k >= len(f4.ROWS):
        ax.axis('off')
        continue
    region = f4.ROWS[k]
    i = 0
    for a, b in SETS:
        for plat, _ in f4.PLATFORMS:
            tot = sum(n[(f'{a}vs{b}', plat, region, c)] for c, _, _ in CLASSES)
            left = 0.0
            for c, _, col in CLASSES:
                w = 100 * n[(f'{a}vs{b}', plat, region, c)] / tot
                ax.barh(ypos[i], w, left=left, height=0.8, color=col, edgecolor='white', lw=0.3)
                left += w
            i += 1
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100])
    ax.set_title(f4.LAB[region], fontweight='bold')
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=5.6)
    ax.tick_params(axis='y', length=0)
axes.flat[-1].legend(handles=[Patch(color=col, label=lab) for _, lab, col in CLASSES], loc='center', fontsize=6.3,
                     title='Allele frequency\nin gnomAD genomes', title_fontsize=6.3)
fig.supxlabel('Share of the pipeline-unique variants of the cohort (%)', fontsize=7, y=0.04)
fig.subplots_adjust(wspace=0.2, hspace=0.35)
f4.supp_title(fig, 'Allele frequency in gnomAD genomes of the pipeline-unique variants, by region')
f4.save(fig, 'SuppFig_uniquevariants_frequency_classes')
print('SUPP_UV_FREQ_CLASSES_DONE')
