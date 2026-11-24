#!/usr/bin/env python3
"""Figure 3 in one file (the panels of figure3_panels.py): PDF, and PNG and TIFF at 600 dpi (TIFF RGB with LZW). Drawn at the printed
size, 183 x 216 mm, and saved without cropping; the rows are placed in inches from the top.
Usage: python3 figure3_combine.py [--c rows|columns|dots]
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure3_panels as f3  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402

f3.style()
plt.rcParams.update({'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'legend.fontsize': 7, 'axes.linewidth': 1.0, 'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
                     'pdf.fonttype': 42})
D, T = f3.load()
cs = f3.c_style()
cols = cs == 'columns'
W, H = 7.205, 8.50                              # 183 mm wide, at most 220 mm tall: the printed size
L, R = 0.165, 0.9                               # a and b: the row labels on the left, the columns of marks of b on the right
LC, RC = (0.10, 0.985) if cols else (L, R)      # c (columns): no marks on the right, so wider axes for the 7-pt marks above it
fig = plt.figure(figsize=(W, H))


def block(top, height, left, right, wspace):
    """three axes sharing the y axis, top and height in inches from the top of the figure"""
    sub = fig.add_gridspec(1, 3, left=left, right=right, top=1 - top / H, bottom=1 - (top + height) / H, wspace=wspace)
    first = fig.add_subplot(sub[0, 0])
    return [first] + [fig.add_subplot(sub[0, j], sharey=first) for j in (1, 2)]


# rows, in inches from the top: tag and title, then the axes, the x axis and the legend
TAG = {'a': 0.04, 'b': 2.95, 'c': 5.62}
blocks = [block(0.44, 1.5, L, R, 0.14), block(3.42, 1.5, L, R, 0.45), block(6.29, 1.5, LC, RC, 0.08 if cols else 0.45)]
LEG = {0: 2.37, 1: 5.36, 2: 8.45}
for ax, (typ, tl) in zip(blocks[0], f3.TYPES):
    f3.draw_pairs(ax, D, typ, ax is blocks[0][0])
    ax.set_title(tl, fontweight='bold')
blocks[0][1].set_xlabel('Genotype concordance in the genome (%)')
for ax, (typ, tl) in zip(blocks[1], f3.TYPES):
    f3.draw_same_pipeline_vs_platform(ax, T, typ, ax is blocks[1][0])
    ax.set_title(tl, fontweight='bold', pad=11)
blocks[1][1].set_xlabel('Same pipeline across platforms minus same platform across pipelines (percentage points)')
for ax, (typ, tl) in zip(blocks[2], f3.TYPES):
    f3.draw_across_platforms(ax, D, T, typ, cs, ax is blocks[2][0])
    ax.set_title(tl, fontweight='bold', pad=24 if cols else 11)
if cols:
    blocks[2][0].set_ylabel('Genotype concordance\nacross platforms (%)')
else:
    blocks[2][1].set_xlabel('Genotype concordance across platforms (%)')
for i, handles, ncol in ((0, f3.legend_pairs(), 3), (1, f3.legend_series(), 2), (2, f3.legend_pipelines(cs != 'dots'), 3)):
    x = (L + R) / 2 if i < 2 else (LC + RC) / 2
    fig.legend(handles=handles, loc='lower center' if i == 2 else 'upper center', bbox_to_anchor=(x, 1 - LEG[i] / H), ncol=ncol,
               handlelength=2.0, columnspacing=1.2, borderaxespad=0, borderpad=0)
for tag in 'abc':
    y = 1 - TAG[tag] / H
    fig.text(0.008, y, tag, fontsize=10, fontweight='bold', va='top')
    fig.text(0.043, y - 0.001, f3.PANEL_TITLE[tag], fontsize=8, fontweight='bold', va='top')
fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print(f'drawn content {bb.x0 * 25.4:.1f}-{bb.x1 * 25.4:.1f} x {(H - bb.y1) * 25.4:.1f}-{(H - bb.y0) * 25.4:.1f} mm of a '
      f'{W * 25.4:.0f} x {H * 25.4:.0f} mm page')
os.makedirs(os.path.join(f3.HERE, 'figures'), exist_ok=True)
for ext, kw in (('pdf', {}), ('png', {'dpi': 600}), ('tiff', {'dpi': 600})):
    fig.savefig(os.path.join(f3.HERE, 'figures', f'Figure3.{ext}'), facecolor='white', **kw)
tif = os.path.join(f3.HERE, 'figures', 'Figure3.tiff')
Image.open(tif).convert('RGB').save(tif, compression='tiff_lzw', dpi=(600, 600))   # flatten the alpha channel
print('written figures/Figure3')
print('FIGURE3_DONE')
