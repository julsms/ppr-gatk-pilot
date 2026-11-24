#!/usr/bin/env python3
"""Figure 7: choice of the VQSR levels for the GATK call sets on HG002 and HG001, Illumina left and Aviti right. a F1 by region at
SNP 99.0, 99.7, 99.9 (INDEL 99.0) and INDEL 99.0, 99.7, 99.9 (SNP 99.9); b HG001 with the cohort and alone at 99.9/99.9; c false
positives and false negatives with minus without the ExcessHet filter, HG001 with the cohort at 99.9/99.9. Drawn at the printed
size, 183 x 218 mm, and saved without cropping.
Input: figure7_table.tsv (figure7_table.py). Output: figures/Figure7.pdf, .png and .tiff (DPI, default 600; TIFF RGB with LZW).
"""
import csv
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402
from PIL import Image  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DPI = int(os.environ.get('DPI', '600'))
INK, GRID, MUTED = '#1f2933', '#e3e7eb', '#7b8794'
COL = {('HG001', 'Illumina'): '#353F56', ('HG001', 'Aviti'): '#7A7977'}
FPC, FNC = '#998829', '#567ABE'
REGIONS = ['BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']
# the level pairs that change one threshold and keep the other: SNP 99.0, 99.7, 99.9 at INDEL 99.0; INDEL 99.0, 99.7, 99.9 at SNP 99.9
LEVELS = {'SNP': [('99.0', 'BP 99.0/99.0'), ('99.7', 'BP 99.7/99.0'), ('99.9', 'BP 99.9/99.0')],
          'INDEL': [('99.0', 'BP 99.9/99.0'), ('99.7', 'BP 99.9/99.7'), ('99.9', 'BP 99.9/99.9')]}
SHADE = '#f4f6f8'
A_ROWS = [('SNP', 'HG002'), ('SNP', 'HG001'), ('INDEL', 'HG002'), ('INDEL', 'HG001')]
PLATS = ('Illumina', 'Aviti')
D_LINE = {'99.9': dict(color='#1A4B8B', marker='o', ms=4.0, ls='-', lw=1.6, mfc='#1A4B8B', z=3),
          '99.7': dict(color='#6A6966', marker='D', ms=5.0, ls=(0, (3.2, 1.8)), lw=1.1, mfc='white', z=4),
          '99.0': dict(color='#C0AA1E', marker='v', ms=6.0, ls=(0, (1.0, 1.4)), lw=1.3, mfc='white', z=5)}

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.0, 'axes.labelsize': 8.0, 'xtick.labelsize': 7.0,
                     'ytick.labelsize': 7.0, 'axes.edgecolor': '#9aa5b1', 'axes.labelcolor': INK, 'xtick.color': INK,
                     'ytick.color': INK, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
                     'axes.unicode_minus': False, 'axes.linewidth': 1.0, 'xtick.major.width': 1.0, 'ytick.major.width': 1.0})

with open(os.path.join(HERE, 'figure7_table.tsv'), newline='') as fh:
    ROWS = list(csv.DictReader(fh, delimiter='\t'))
IDX = {(r['genome'], r['platform'], r['calling'], r['setting'], r['region'], r['type'], r['filter']): r for r in ROWS}


def get(g, p, c, s, region, t, filt='PASS'):
    return IDX.get((g, p, c, s, region, t, filt))


def f1(g, p, c, set_, region, t):
    r = get(g, p, c, set_, region, t)
    return float(r['F1']) if r else None


def levels_in(g, p, t):
    """the thresholds benchmarked for this genome, platform and type (calls joint with the cohort), in the order 99.0, 99.7, 99.9"""
    return [(lev, set_) for lev, set_ in LEVELS[t] if get(g, p, 'cohort', set_, 'BB', t)]


def region_axis(ax, labels):
    ax.axvspan(-0.5, 0.5, color=SHADE, lw=0, zorder=0)
    ax.set_xlim(-0.5, len(REGIONS) - 0.5)
    ax.set_xticks(range(len(REGIONS)))
    if labels:
        ax.set_xticklabels(REGIONS, rotation=45, ha='right', rotation_mode='anchor', fontsize=7)
    else:
        ax.set_xticklabels([])
    ax.tick_params(axis='x', length=0)
    ax.grid(axis='x', color=GRID, lw=0.5, zorder=0)


class Inch:
    """axes and text placed in inches from the left and from the top of a W x H inch figure"""
    def __init__(self, fig, W, H):
        self.fig, self.W, self.H = fig, W, H

    def ax(self, x, y, w, h):
        return self.fig.add_axes([x / self.W, (self.H - y - h) / self.H, w / self.W, h / self.H])

    def text(self, x, y, s, **kw):
        return self.fig.text(x / self.W, (self.H - y) / self.H, s, **kw)

    def at(self, x, y):
        return x / self.W, (self.H - y) / self.H


def f1_floor():
    """per variant type: the lower end of the F1 axis of panel a (a multiple of 0.05 below the lowest F1 of the four rows of that type)"""
    lo = {}
    for t in ('SNP', 'INDEL'):
        vals = []
        for _, g in [r for r in A_ROWS if r[0] == t]:
            for p in PLATS:
                lv = dict(levels_in(g, p, t))
                if '99.9' not in lv:
                    continue
                for set_ in lv.values():
                    vals.extend(f1(g, p, 'cohort', set_, region, t) for region in REGIONS)
        lo[t] = np.floor((min(vals) - 0.02) * 20) / 20
    return lo


def header(L, x0, colw, gap, y, size=8):
    for j, p in enumerate(PLATS):
        L.text(x0 + j * (colw + gap) + colw / 2, y, p, fontsize=size, fontweight='bold', ha='center', va='bottom')


def panel_a(L, x0, y0, colw, gap, ph, rowgap, blockgap):
    """the line panels of a; returns the bottom of the last row (inches from the top)"""
    X = np.arange(len(REGIONS))
    floor = f1_floor()
    y = y0
    for k, (t, g) in enumerate(A_ROWS):
        lo = floor[t]
        if k == 2:
            header(L, x0, colw, gap, y - 0.04)
        for j, p in enumerate(PLATS):
            ax = L.ax(x0 + j * (colw + gap), y, colw, ph)
            region_axis(ax, k in (1, 3))
            ax.grid(axis='y', color='#cfd6dc', lw=0.6, zorder=0)
            ax.set_ylim(lo, 1.01)
            ax.set_yticks(np.arange(np.ceil(lo * 5 - 1e-9) / 5, 1.0001, 0.2))
            ax.yaxis.set_major_formatter(lambda v, _: f'{v:.1f}')
            if j == 0:
                ax.set_ylabel(f'{g}\n{"SNV" if t == "SNP" else "indel"} F1')
            else:
                ax.tick_params(labelleft=False)
            for lev, set_ in sorted(levels_in(g, p, t), key=lambda z: z[0] != '99.9'):
                st = D_LINE[lev]
                ys = [f1(g, p, 'cohort', set_, region, t) for region in REGIONS]
                ax.plot(X, ys, color=st['color'], ls=st['ls'], lw=st['lw'], marker=st['marker'], ms=st['ms'], mfc=st['mfc'],
                        mec=st['color'], mew=1.0, zorder=st['z'], clip_on=False)
        y += ph + (blockgap if k in (1, 3) else rowgap)
    return y - blockgap + 0.5


def panel_b(L, x0, y0, colw, gap, ph):
    """HG001, SNVs and indels together at 99.9/99.9, with the cohort (solid) and alone (dashed), a platform in each column"""
    X = np.arange(len(REGIONS))
    vals = [f1('HG001', p, c, 'BP 99.9/99.9', r, 'ALL') for p in PLATS for c in ('cohort', 'single') for r in REGIONS]
    lo = np.floor((min(vals) - 0.01) * 20) / 20
    for j, p in enumerate(PLATS):
        ax = L.ax(x0 + j * (colw + gap), y0, colw, ph)
        region_axis(ax, True)
        ax.grid(axis='y', color='#cfd6dc', lw=0.6, zorder=0)
        ax.set_ylim(lo, 1.005)
        ax.set_yticks(np.arange(np.ceil(lo * 10 - 1e-9) / 10, 1.0001, 0.1))
        ax.yaxis.set_major_formatter(lambda v, _: f'{v:.1f}')
        col = COL[('HG001', p)]
        for c, ls, mfc, lw, z in (('single', (0, (3.2, 1.8)), 'white', 1.1, 4), ('cohort', '-', col, 1.6, 3)):
            ax.plot(X, [f1('HG001', p, c, 'BP 99.9/99.9', r, 'ALL') for r in REGIONS], color=col, ls=ls, lw=lw, marker='s',
                    ms=4.0 if c == 'cohort' else 4.8, mfc=mfc, mec=col, mew=1.0, zorder=z, clip_on=False)
        if j == 0:
            ax.set_ylabel('HG001\nSNV and\nindel F1')
            hand = [Line2D([], [], color=MUTED, lw=1.6, marker='s', ms=4.0, label='With the cohort'),
                    Line2D([], [], color=MUTED, lw=1.1, ls=(0, (3.2, 1.8)), marker='s', ms=4.8, mfc='white', label='Alone')]
            ax.legend(handles=hand, loc='lower left', fontsize=7, handlelength=2.6, handletextpad=0.4, borderaxespad=0.4,
                      labelspacing=0.3)
        else:
            ax.tick_params(labelleft=False)
    return y0 + ph + 0.55


def panel_c(L, x0, y0, colw, gap, ph, rowgap):
    """ExcessHet filter: with minus without, false positives and false negatives by region; SNVs (top row) and indels (bottom row)"""
    X = np.arange(len(REGIONS))
    data = {}
    for t in ('SNP', 'INDEL'):
        for p in PLATS:
            fp, fn = [], []
            for region in REGIONS:
                a = get('HG001', p, 'cohort', 'BP 99.9/99.9 no ExcessHet', region, t)
                b = get('HG001', p, 'cohort', 'BP 99.9/99.9', region, t)
                fp.append(int(b['FP']) - int(a['FP']))
                fn.append(int(b['FN']) - int(a['FN']))
            data[(t, p)] = (np.array(fp), np.array(fn))
    y = y0
    left = []
    for r, t in enumerate(('SNP', 'INDEL')):
        every = [v for p in PLATS for v in (*data[(t, p)], data[(t, p)][0] + data[(t, p)][1])]
        lo, hi = min(0, min(v.min() for v in every)), max(0, max(v.max() for v in every))
        pad = 0.08 * (hi - lo)
        for j, p in enumerate(PLATS):
            ax = L.ax(x0 + j * (colw + gap), y, colw, ph)
            fp, fn = data[(t, p)]
            ax.bar(X - 0.18, fp, width=0.34, color=FPC, lw=0, zorder=3)
            ax.bar(X + 0.18, fn, width=0.34, color=FNC, lw=0, zorder=3)
            ax.plot(X, fp + fn, marker='D', color=INK, ms=4.0, ls='', zorder=4)
            ax.axhline(0, color='#52606d', lw=1.0, zorder=4)
            region_axis(ax, r == 1)
            ax.grid(axis='y', color='#cfd6dc', lw=0.6, zorder=0)
            ax.set_ylim(lo - pad, hi + pad)
            ax.yaxis.set_major_locator(MaxNLocator(4, steps=[1, 2, 5, 10]))
            ax.yaxis.set_major_formatter(lambda v, _: f'{v:,.0f}')
            if j == 0:
                ax.set_ylabel('Change in\n' + ('SNV' if t == 'SNP' else 'indel') + ' errors')
                left.append(ax)
            else:
                ax.tick_params(labelleft=False)
        y += ph + (0.14 if r == 0 else rowgap)
    for ax in left:          # the two y titles in one column although the tick labels differ in width
        ax.yaxis.set_label_coords(-0.155, 0.5)
    return y - 0.14 + 0.5


def figure7():
    W, H = 183 / 25.4, 8.6          # the printed size, saved without cropping
    fig = plt.figure(figsize=(W, H))
    L = Inch(fig, W, H)
    x0, colw, gap = 0.84, 3.09, 0.12
    L.text(0.04, 0.08, 'a', fontsize=10, fontweight='bold', va='top')
    L.text(0.24, 0.09, 'F1 by region at the tested VQSR levels (joint-called with the cohort)', fontsize=8, fontweight='bold', va='top')
    h = [Line2D([], [], color=D_LINE[lev]['color'], ls=D_LINE[lev]['ls'], lw=D_LINE[lev]['lw'], marker=D_LINE[lev]['marker'],
                ms=D_LINE[lev]['ms'], mfc=D_LINE[lev]['mfc'], mec=D_LINE[lev]['color'], mew=1.0, label=lab)
         for lev, lab in (('99.9', '99.9 (chosen)'), ('99.7', '99.7'), ('99.0', '99.0'))]
    fig.legend(handles=h, loc='upper left', bbox_to_anchor=L.at(0.19, 0.28), ncol=3, fontsize=7, handlelength=3.0, columnspacing=1.6,
               borderaxespad=0, title='VQSR truth sensitivity level that changes', title_fontsize=7, alignment='left')
    L.text(3.0, 0.32, 'Top two rows, SNVs: the SNP level changes, INDEL is fixed at 99.0.\n'
           'Bottom two rows, indels: the INDEL level changes, SNP is fixed at 99.9.', fontsize=7, ha='left', va='top', color=INK,
           linespacing=1.4)
    ya = 0.82
    header(L, x0, colw, gap, ya - 0.04)
    bottom = panel_a(L, x0, ya, colw, gap, 0.60, 0.12, 0.55)
    yb = bottom + 0.01
    L.text(0.04, yb, 'b', fontsize=10, fontweight='bold', va='top')
    L.text(0.24, yb + 0.01, 'HG001 alone or with the cohort (SNVs and indels together, 99.9/99.9)', fontsize=8, fontweight='bold', va='top')
    header(L, x0, colw, gap, yb + 0.30)
    yc = panel_b(L, x0, yb + 0.34, colw, gap, 0.72)
    L.text(0.04, yc, 'c', fontsize=10, fontweight='bold', va='top')
    L.text(0.24, yc + 0.01, 'ExcessHet filter in HG001 with the cohort: calls with the filter minus calls without it', fontsize=8,
           fontweight='bold', va='top')
    hc = [Patch(facecolor=FPC, label='False positives'), Patch(facecolor=FNC, label='False negatives'),
          Line2D([], [], marker='D', color=INK, ls='', ms=4.0, label='Net change')]
    fig.legend(handles=hc, loc='upper left', bbox_to_anchor=L.at(0.19, yc + 0.20), ncol=3, fontsize=7, borderaxespad=0,
               columnspacing=1.6, handlelength=1.4)
    header(L, x0, colw, gap, yc + 0.46)
    panel_c(L, x0, yc + 0.50, colw, gap, 0.66, 0.12)
    bad = [x.get_text() for x in fig.findobj(plt.Text) if any(ord(ch) in (0x2012, 0x2013, 0x2014, 0x2015, 0x2212) for ch in x.get_text())]
    if bad:
        sys.exit(f'figure7: long dash in labels: {bad[:3]}')
    return fig


fig = figure7()
os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True)
for ext, kw in (('pdf', {}), ('png', {'dpi': DPI}), ('tiff', {'dpi': DPI, 'pil_kwargs': {'compression': 'tiff_lzw'}})):
    fig.savefig(os.path.join(HERE, 'figures', f'Figure7.{ext}'), facecolor='white', **kw)
plt.close(fig)
tif = os.path.join(HERE, 'figures', 'Figure7.tiff')     # RGB, no alpha channel
Image.open(tif).convert('RGB').save(tif, dpi=(DPI, DPI), compression='tiff_lzw')
print(f'written figures/Figure7 (PNG and TIFF at {DPI} dpi)')
print('FIGURE7_DONE')
