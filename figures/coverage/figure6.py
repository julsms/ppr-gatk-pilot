#!/usr/bin/env python3
"""figure6.py - Figure 6, drawn and saved at the printed size, 183 x 165 mm: a depth on Illumina and Aviti; b concordance within a
platform against depth; c concordance between platforms against the lower of the two depths; d k-mer support against depth.
Output: figures/Figure6.{pdf,png,tiff} (PNG and TIFF at DPI dots per inch, default 600; TIFF RGB with LZW).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure6_panels as f6  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

f6.style()
plt.rcParams.update({'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'axes.linewidth': 1.0,
                     'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 'pdf.fonttype': 42})
W, H = 7.205, 6.5                  # 183 mm wide: the printed size
fig = plt.figure(figsize=(W, H))
t1, b1, t2, b2 = 0.935, 0.615, 0.445, 0.075
ax_a = fig.add_subplot(fig.add_gridspec(1, 1, left=0.075, right=0.30, top=t1, bottom=b1)[0])
gb = fig.add_gridspec(1, 2, left=0.40, right=0.985, top=t1, bottom=b1, wspace=0.08)
ax_b = [fig.add_subplot(gb[0])]
ax_b.append(fig.add_subplot(gb[1], sharey=ax_b[0]))
ax_c = fig.add_subplot(fig.add_gridspec(1, 1, left=0.075, right=0.36, top=t2, bottom=b2)[0])
gd = fig.add_gridspec(1, 2, left=0.46, right=0.985, top=t2, bottom=b2, wspace=0.25)
ax_d = [fig.add_subplot(gd[0]), fig.add_subplot(gd[1])]
f6.draw_a(ax_a)
f6.draw_b(ax_b)
ax_b[1].tick_params(axis='y', labelleft=False)
f6.draw_c(ax_c)
f6.draw_d(ax_d)
for x, y, tag in ((0.005, 0.995, 'a'), (0.33, 0.995, 'b'), (0.005, 0.525, 'c'), (0.39, 0.525, 'd')):
    fig.text(x, y, tag, fontsize=10, fontweight='bold', va='top', ha='left')
    fig.text(x + 0.03, y - 0.003, f6.TITLES[tag], fontsize=8, fontweight='bold', va='top', ha='left')
f6.save(fig, 'Figure6', tiff=True, tight=False)
w, h = W * 25.4, H * 25.4
fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print(f'drawn content {bb.x0 * 25.4:.1f}-{bb.x1 * 25.4:.1f} x {(H - bb.y1) * 25.4:.1f}-{(H - bb.y0) * 25.4:.1f} mm')
plt.close(fig)
print(f'written figures/Figure6 ({w:.0f} x {h:.0f} mm; PNG and TIFF at {f6.DPI} dpi)')
print('FIGURE6_DONE')
