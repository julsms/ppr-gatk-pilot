#!/usr/bin/env python3
"""Figure 1: panel a (SNVs, left) and panel b (indels, right), drawn with load() and draw_column() of figure1_panels.py, one legend,
panel tags with a short title. Drawn at the printed size, 183 mm wide; the script stops if the content runs outside that width or
the figure is taller than 220 mm.
Inputs: benchmark_data.tsv. Outputs: figures/Figure1.{pdf,png,tiff} (PNG and TIFF 600 dpi, TIFF RGB with LZW).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure1_panels as f1  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.transforms import Bbox  # noqa: E402
from PIL import Image  # noqa: E402

WIDTH_IN, MAX_HEIGHT_MM = 183 / 25.4, 220
f1.style()
data = f1.load()
fig, axes = plt.subplots(3, 2, figsize=(WIDTH_IN, 7.6), sharex=True)
f1.draw_column(axes[:, 0], data, 'SNP', show_ylabel=True)
f1.draw_column(axes[:, 1], data, 'INDEL', show_ylabel=False)
TITLE = {'a': 'Precision, recall and F1 of SNVs', 'b': 'Precision, recall and F1 of indels'}
for ax, tag in ((axes[0, 0], 'a'), (axes[0, 1], 'b')):
    ax.set_title('')   # the column title (SNVs, Indels) is in the panel title
    ax.text(-0.14, 1.06, tag, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')
    ax.text(-0.07, 1.065, TITLE[tag], transform=ax.transAxes, fontsize=8, fontweight='bold', va='bottom')
fig.legend(handles=f1.legend_handles(), loc='upper center', ncol=4, bbox_to_anchor=(0.53, 0.035), fontsize=7,
           columnspacing=1.2, handlelength=2.6)
fig.subplots_adjust(left=0.095, right=0.99, top=0.955, bottom=0.11, hspace=0.12, wspace=0.17)
# saved 183 mm wide whatever the text needs: the crop keeps the full width and trims only the height
fig.canvas.draw()
tb = fig.get_tightbbox(fig.canvas.get_renderer())
if tb.x0 < 0 or tb.x1 > WIDTH_IN:
    raise SystemExit(f'Figure 1 content runs outside the 183 mm width: {tb.x0:.3f} to {tb.x1:.3f} in')
box = Bbox([[0, tb.y0 - 0.03], [WIDTH_IN, tb.y1 + 0.03]])
if box.height * 25.4 > MAX_HEIGHT_MM:
    raise SystemExit(f'Figure 1 is {box.height * 25.4:.0f} mm tall (at most {MAX_HEIGHT_MM})')
out = os.path.join(f1.HERE, 'figures', 'Figure1')
os.makedirs(os.path.dirname(out), exist_ok=True)
for ext, kw in (('pdf', {}), ('png', {'dpi': 600})):
    fig.savefig(f'{out}.{ext}', bbox_inches=box, facecolor='white', **kw)
Image.open(f'{out}.png').convert('RGB').save(f'{out}.tiff', compression='tiff_lzw', dpi=(600, 600))
print(f'written figures/Figure1.pdf .png .tiff ({box.width * 25.4:.0f} x {box.height * 25.4:.0f} mm)')
print('FIGURE1_DONE')
