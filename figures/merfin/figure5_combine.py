#!/usr/bin/env python3
"""figure5_combine.py - Figure 5, drawn and saved at the printed size, 183 x 214 mm: a k-mer support in the cohort; b HG002, k-mer
support and GIAB precision; c HG002, k-mer context of the supported records; d HG002 Illumina, k-mer support and GIAB precision by
region, with the differences in one row above the data.
Output: figures/Figure5.{pdf,png,tiff} (PNG and TIFF at DPI dots per inch, default 600; TIFF RGB with LZW).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure5_panels as f5  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402

DPI = int(os.environ.get('DPI', '600'))
f5.style()
plt.rcParams.update({'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'legend.fontsize': 7, 'axes.linewidth': 1.0, 'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
                     'pdf.fonttype': 42})
R, NV, HG, HR = f5.load()
W, H = 7.205, 8.42                  # 183 mm wide, at most 220 mm tall: the printed size
fig = plt.figure(figsize=(W, H))
TITLE = 0.44 / H                     # room for the tag and title above a row, over the titles of the axes
GAP = 0.12 / H                       # between a legend and the next row


def lowest(artists):
    """the lowest point (figure fraction) of some axes with their tick labels"""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    return min(a.get_tightbbox(rend).y0 for a in artists) / fig.bbox.height


def legend_under(artists, handles, ncol, xc):
    """a legend just under the artists; returns the bottom of the legend (figure fraction)"""
    leg = fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(xc, lowest(artists) - 0.03 / H), ncol=ncol, fontsize=7,
                     borderaxespad=0, handletextpad=0.5, columnspacing=1.4)
    fig.canvas.draw()
    return leg.get_window_extent(fig.canvas.get_renderer()).y0 / fig.bbox.height


# a (left) and b (right)
t1 = 1 - TITLE
b1 = t1 - 1.8 / H
ax_a = fig.add_subplot(fig.add_gridspec(1, 1, left=0.2, right=0.36, top=t1, bottom=b1 - 0.4 / H)[0])
g1 = fig.add_gridspec(1, 2, left=0.45, right=0.995, top=t1, bottom=b1, wspace=0.05)
ax_b = [fig.add_subplot(g1[0])]
ax_b.append(fig.add_subplot(g1[1], sharey=ax_b[0]))
ha = f5.draw_support_bars(ax_a, R, fs=7)
hb = f5.draw_merfin_vs_giab(ax_b, HG, fs=7, vertical=True, col=f5.OLIVE)
ax_b[1].tick_params(axis='y', labelleft=False)
ax_b[0].set_ylabel('Records inside the\nGIAB benchmark (%)', fontsize=8)
y = legend_under(ax_b, hb, 2, 0.72)
y = min(y, legend_under([ax_a], ha, 2, 0.20))
# c: back-to-back bars; the set names on the left axis need a wider margin than the other rows
t2 = y - GAP - TITLE
g2 = fig.add_gridspec(1, 2, left=0.2, right=0.99, top=t2, bottom=t2 - 1.45 / H, wspace=0.08)
ax_c = [fig.add_subplot(g2[0])]
ax_c.append(fig.add_subplot(g2[1], sharey=ax_c[0]))
hc = f5.draw_context(ax_c, HG, fs=7)
ax_c[1].tick_params(axis='y', labelleft=False)
y = legend_under(ax_c, hc, 3, 0.54)
# d
t3 = y - GAP - TITLE
g3 = fig.add_gridspec(1, 2, left=0.09, right=0.995, top=t3, bottom=t3 - 1.25 / H, wspace=0.06)
ax_d = [fig.add_subplot(g3[0])]
ax_d.append(fig.add_subplot(g3[1], sharey=ax_d[0]))
hd = f5.draw_regions_lines(ax_d, HR, 'illumina', side_by_side=True, fs=7, col=f5.OLIVE)
ax_d[1].tick_params(axis='y', labelleft=False)
ax_d[0].set_ylabel('Records inside the\nGIAB benchmark (%)', fontsize=8)
legend_under(ax_d, hd, 2, 0.54)
for x, top, tag in ((0.005, t1, 'a'), (0.375, t1, 'b'), (0.005, t2, 'c'), (0.005, t3, 'd')):
    fig.text(x, top + TITLE * 0.95, tag, fontsize=10, fontweight='bold', va='top', ha='left')
    fig.text(x + 0.03, top + TITLE * 0.95 - 0.004, f5.PANEL_TITLE[tag], fontsize=8, fontweight='bold', va='top', ha='left')
os.makedirs(os.path.join(f5.HERE, 'figures'), exist_ok=True)
fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print(f'drawn content {bb.x0 * 25.4:.1f}-{bb.x1 * 25.4:.1f} x {(H - bb.y1) * 25.4:.1f}-{(H - bb.y0) * 25.4:.1f} mm of a '
      f'{W * 25.4:.0f} x {H * 25.4:.0f} mm page')
for ext, kw in (('pdf', {}), ('png', {'dpi': DPI}), ('tiff', {'dpi': DPI})):
    fig.savefig(os.path.join(f5.HERE, 'figures', f'Figure5.{ext}'), facecolor='white', **kw)
tif = os.path.join(f5.HERE, 'figures', 'Figure5.tiff')
Image.open(tif).convert('RGB').save(tif, compression='tiff_lzw', dpi=(DPI, DPI))   # flatten the alpha channel
print(f'written figures/Figure5 ({W * 25.4:.0f} x {H * 25.4:.0f} mm, {DPI} dpi)')
print('FIGURE5_DONE')
