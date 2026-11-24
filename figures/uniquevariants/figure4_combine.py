#!/usr/bin/env python3
"""figure4_combine.py - Figure 4 (variants found by one pipeline and not by the other), panels a to d, drawn with the functions of
figure4_panels.py, at the printed size, 183 x 201 mm, and saved without cropping. Outputs: figures/Figure4.{pdf,png,tiff} (PNG and
TIFF at 600 dpi, TIFF RGB with LZW)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure4_panels as f4  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402

f4.style()
# the text and line sizes of the main figures at 183 mm; the supplementary figures keep the sizes of f4.style()
plt.rcParams.update({'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7, 'axes.linewidth': 1.0,
                     'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
C, H, _ = f4.load()
_, N = f4.load_r27()
WIDTH, HEIGHT = 183 / 25.4, 7.9          # the printed size
fig = plt.figure(figsize=(WIDTH, HEIGHT))
# rows: a, its legend, b, c, its legend, d, its legend; b to d keep an empty top part for their titles
rows = fig.add_gridspec(7, 1, height_ratios=[2.4, 0.25, 1.6, 1.6, 0.22, 1.0, 0.2], hspace=0.38, left=0.12,
                        right=0.98, top=1 - 0.52 / HEIGHT, bottom=0.29 / HEIGHT)
# a: two rows of violins (GATK vs PPR, DV vs PPR), three groups of regions each; region names under the lower row only
a_rows = rows[0].subgridspec(2, 3, width_ratios=[6, 3, 2], wspace=0.22, hspace=0.22)
a_axes = []
for i, pair in enumerate(f4.PAIRS):
    axes = [fig.add_subplot(a_rows[i, j]) for j in range(3)]
    f4.draw_violins(axes, C, pair, labels=(i == len(f4.PAIRS) - 1))
    axes[0].text(0.0, 1.02, f'{pair[0]} vs {pair[1]}', transform=axes[0].transAxes, fontweight='bold', fontsize=7, va='bottom')
    a_axes.append(axes)
a_axes[0][0].set_ylabel('Pipeline-unique variants\nper individual', y=-0.11)
lax = fig.add_subplot(rows[1])
lax.axis('off')
lax.legend(handles=f4.legend_violins(), loc='center', ncol=6, fontsize=7, handlelength=1.4, columnspacing=1.0)
# b: gnomAD frequency, four histograms
b_row = rows[2].subgridspec(2, 4, height_ratios=[0.25, 1.2], hspace=0, wspace=0.36)
b_axes = []
for j, (pair, (plat, pl)) in enumerate([(p, q) for p in f4.PAIRS for q in f4.PLATFORMS]):
    ax = fig.add_subplot(b_row[1, j])
    f4.draw_af(ax, H, pair, plat, ylabel=(j == 0))
    ax.set_title(f'{pair[0]} vs {pair[1]}, {pl}', fontweight='bold', fontsize=7)
    ax.set_xlabel('AF in gnomAD genomes', fontsize=8, labelpad=1)
    b_axes.append(ax)
b_axes[0].set_ylabel('Pipeline-unique\nvariants found\nin gnomAD')      # the words of draw_af on three lines, to fit at 8 pt
# c: novel alleles, regions in the three groups of a (lollipops, own axes); d: the caller step and the graph step by region
c_spec = rows[3].subgridspec(2, 3, height_ratios=[0.5, 1.1], width_ratios=[6, 3, 2], hspace=0, wspace=0.22)
c_axes = [fig.add_subplot(c_spec[1, j]) for j in range(3)]
f4.draw_novel(c_axes, N)
d_ax = fig.add_subplot(rows[5].subgridspec(2, 1, height_ratios=[0.2, 0.75], hspace=0)[1, 0])
f4.draw_caller_graph(d_ax, C, 'ALL')
fig.align_ylabels([c_axes[0], d_ax])
for row, handles in ((rows[4], f4.legend_novel()), (rows[6], f4.legend_caller_graph())):
    lax = fig.add_subplot(row)
    lax.axis('off')
    lax.legend(handles=handles, loc='center', ncol=len(handles), fontsize=7, handlelength=1.4, columnspacing=0.9)
# tag offsets of a 9.7-inch layout, scaled to keep the same distance in inches
for ax, tag in ((a_axes[0][0], 'a'), (b_axes[0], 'b'), (c_axes[0], 'c'), (d_ax, 'd')):
    y = ax.get_position().y1 + {'a': 0.026, 'b': 0.032}.get(tag, 0.014) * 9.7 / HEIGHT
    fig.text(0.01, y, tag, fontsize=10, fontweight='bold', va='bottom')
    fig.text(0.04, y + 0.001, f4.PANEL_TITLE[tag], fontsize=8, fontweight='bold', va='bottom')
os.makedirs(os.path.join(f4.HERE, 'figures'), exist_ok=True)
for ext, kw in (('pdf', {}), ('png', {'dpi': 600}), ('tiff', {'dpi': 600, 'pil_kwargs': {'compression': 'tiff_lzw'}})):
    fig.savefig(os.path.join(f4.HERE, 'figures', f'Figure4.{ext}'), facecolor='white', **kw)
tif = os.path.join(f4.HERE, 'figures', 'Figure4.tiff')     # RGB, no alpha channel
Image.open(tif).convert('RGB').save(tif, dpi=(600, 600), compression='tiff_lzw')
print('written figures/Figure4')
print('FIGURE4_DONE')
