#!/usr/bin/env python3
"""Figure 2: panels a to e drawn with the functions of figure2_panels.py, the legends at the bottom, panel tags with a short title.
Drawn at the printed size, 183 mm wide; the script stops if the content runs outside that width or the figure is taller than 220 mm.
Inputs: graph_arms_data.tsv, graph_arms_roc.tsv, graph_arms_mapping.tsv. Outputs: figures/Figure2.{pdf,png,tiff} (PNG and TIFF 600 dpi,
TIFF RGB with LZW).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure2_panels as f2  # noqa: E402
import graph_setups as gs  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.transforms import Bbox  # noqa: E402
from PIL import Image  # noqa: E402

WIDTH_IN, MAX_HEIGHT_MM = 183 / 25.4, 220
f2.style()
data = f2.load()
fig = plt.figure(figsize=(WIDTH_IN, 7.0))
# left margin wide enough for the set-up names of panel e; panel e has one row per set-up (five)
rows = fig.add_gridspec(3, 1, height_ratios=[2.3, 2.3, 1.3], hspace=0.62, left=0.125, right=0.985, top=0.96, bottom=0.04)
top, mid = rows[0].subgridspec(1, 2, wspace=0.25), rows[1].subgridspec(1, 2, wspace=0.25)
low = rows[2].subgridspec(1, len(f2.FIG2_METRICS), wspace=0.1)
ax_a, ax_b = fig.add_subplot(top[0]), fig.add_subplot(top[1])
ax_c, ax_d = fig.add_subplot(mid[0]), fig.add_subplot(mid[1])
ax_e = [fig.add_subplot(low[0, j]) for j in range(len(f2.FIG2_METRICS))]
f2.draw_pr(ax_a, data, 'SNP')
f2.draw_pr(ax_b, data, 'INDEL', show_ylabel=False)
f2.draw_errors(ax_c, data, 'SNP')
f2.draw_errors(ax_d, data, 'INDEL', show_ylabel=False)
f2.draw_alignment(ax_e, data, titles=f2.FIG2_TITLES)
# panel tags in one column on the left (a, c, e) and one in the middle (b, d), each with its title; the column titles (SNVs, Indels)
# of a to d are in the panel titles
TITLE = {'a': 'Precision and recall, SNVs', 'b': 'Precision and recall, indels', 'c': 'Errors, SNVs', 'd': 'Errors, indels',
         'e': 'Alignment metrics, BAM files on chromosomes 1-22, X and Y'}
for ax in (ax_a, ax_b, ax_c, ax_d):
    ax.set_title('')
x_left = ax_a.get_position().x0 - 0.075
x_mid = ax_b.get_position().x0 - 0.06
for ax, x, tag in ((ax_a, x_left, 'a'), (ax_b, x_mid, 'b'), (ax_c, x_left, 'c'), (ax_d, x_mid, 'd'), (ax_e[0], x_left, 'e')):
    y = ax.get_position().y1 + (0.045 if tag == 'e' else 0.016)
    fig.text(x, y, tag, fontsize=10, fontweight='bold', va='bottom')
    fig.text(x + 0.028, y + 0.001, TITLE[tag], fontsize=8, fontweight='bold', va='bottom')
# the legends right under panel e (below its tick labels), the keys right under the set-ups
fig.canvas.draw()
inv = fig.transFigure.inverted()
renderer = fig.canvas.get_renderer()
y_e = min(ax.get_tightbbox(renderer).transformed(inv).y0 for ax in ax_e)
leg = fig.legend(handles=f2.setup_handles(), loc='upper left', bbox_to_anchor=(0.02, y_e - 0.01), fontsize=7, title=gs.LEGEND_TITLE,
                 title_fontsize=7, alignment='left', handletextpad=0.3, labelspacing=0.4, borderaxespad=0)
fig.canvas.draw()
y_leg = leg.get_window_extent(renderer).transformed(inv).y0
fig.legend(handles=f2.key_handles(), loc='upper left', bbox_to_anchor=(0.02, y_leg - 0.006), fontsize=7, ncol=2, alignment='left',
           handlelength=2.6, columnspacing=1.2, labelspacing=0.4, borderaxespad=0)
# saved 183 mm wide whatever the text needs: the crop keeps the full width and trims only the height
fig.canvas.draw()
tb = fig.get_tightbbox(fig.canvas.get_renderer())
if tb.x0 < 0 or tb.x1 > WIDTH_IN:
    raise SystemExit(f'Figure 2 content runs outside the 183 mm width: {tb.x0:.3f} to {tb.x1:.3f} in')
box = Bbox([[0, tb.y0 - 0.03], [WIDTH_IN, tb.y1 + 0.03]])
if box.height * 25.4 > MAX_HEIGHT_MM:
    raise SystemExit(f'Figure 2 is {box.height * 25.4:.0f} mm tall (at most {MAX_HEIGHT_MM})')
out = os.path.join(f2.HERE, 'figures', 'Figure2')
os.makedirs(os.path.dirname(out), exist_ok=True)
for ext, kw in (('pdf', {}), ('png', {'dpi': 600})):
    fig.savefig(f'{out}.{ext}', bbox_inches=box, facecolor='white', **kw)
Image.open(f'{out}.png').convert('RGB').save(f'{out}.tiff', compression='tiff_lzw', dpi=(600, 600))
print(f'written figures/Figure2.pdf .png .tiff ({box.width * 25.4:.0f} x {box.height * 25.4:.0f} mm)')
print('FIGURE2_DONE')
