#!/usr/bin/env python3
"""supp_uniquevariants_titv_hethom.py - Supplementary Figure S18: Ti/Tv of the SNVs and het/hom by region, pooled over the 23
individuals, Illumina and Aviti, for the pairs GATK vs PPR (circles) and DV vs PPR (squares). a: the variants found by one pipeline of
the pair; b: the shared variants plus the variants found by one pipeline of the pair (shared + GATK-only, shared + PPR-only, shared +
DV-only). The shared variants are a grey line in both. Drawn with figure4_panels.draw_titv_hethom.
Method: ratios from the summed counts (ti / tv of the SNVs, het / hom of all records) of the classes shared and unique of
uv_titv_hethom.tsv; shared = the mean of the two sides of a pair, as in figure4_panels.load_r27.
Input: uv_titv_hethom.tsv. Outputs: figures/SuppFig_uniquevariants_titv_hethom.{png,pdf}
"""
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure4_panels as f4  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

PAIRS = ('GATK_PPR', 'DV_PPR')
counts = defaultdict(lambda: np.zeros(4))
for r in f4.read('uv_titv_hethom.tsv'):
    if r['pair'] in PAIRS and r['class'] in ('shared', 'unique') and r['region'] in f4.ROWS:
        counts[(r['pair'], r['set'], r['platform'], r['region'], r['class'])] += [float(r[c]) for c in ('ti', 'tv', 'het', 'hom')]


def ratios(c):
    return (c[0] / c[1] if c[1] else float('nan'), c[2] / c[3] if c[3] else float('nan'))


def table(classes):
    """(pair, set or 'shared', platform, region) -> (Ti/Tv, het/hom) of the summed classes; shared: mean of the two sides"""
    T, shared = {}, defaultdict(list)
    for (pair, s, plat, reg, cl) in list(counts):
        T[(pair, s, plat, reg)] = ratios(sum(counts[(pair, s, plat, reg, c)] for c in classes))
        shared[(pair, 'shared', plat, reg)].append(ratios(counts[(pair, s, plat, reg, 'shared')]))
    for k, v in shared.items():
        T[k] = tuple(float(np.mean([x[i] for x in v])) for i in (0, 1))
    return T


f4.style()
fig = plt.figure(figsize=(7.1, 6.2))
gs = fig.add_gridspec(5, 1, height_ratios=[1, 1, 0.38, 1, 1], hspace=0.12)
items = []
for top, tag, classes, title in ((0, 'a', ('unique',), 'Variants found by one pipeline of the pair'),
                                 (3, 'b', ('shared', 'unique'), 'Shared variants plus the variants found by one pipeline')):
    axes = [fig.add_subplot(gs[top]), fig.add_subplot(gs[top + 1])]
    f4.draw_titv_hethom(axes, table(classes))
    axes[0].tick_params(labelbottom=False)
    items.append((axes[0], tag, title))
fig.legend(handles=f4.legend_r27(), loc='upper center', bbox_to_anchor=(0.5, 0.04), ncol=7, fontsize=5.8)
f4.panel_tags(fig, items)
f4.supp_title(fig, 'Ti/Tv and het/hom by region')
f4.save(fig, 'SuppFig_uniquevariants_titv_hethom')
print('SUPP_UV_TITV_DONE')
