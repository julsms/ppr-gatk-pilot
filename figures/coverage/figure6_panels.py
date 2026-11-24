#!/usr/bin/env python3
"""figure6_panels.py - style, data and drawing functions of Figure 6 and of the depth-by-region supplementary figure (imported), with
the tables of coverage_table.py from this folder. Colours of Figures 1, 3 and 4 (Paul Tol's muted scheme), one marker shape per
series; Spearman rho with Holm-adjusted p in the keys; Theil-Sen guide lines.
"""
import csv
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.text import Text  # noqa: E402
from PIL import Image  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DPI = int(os.environ.get('DPI', '600'))
PLATS = [('illumina', 'Illumina'), ('aviti', 'Aviti')]
INK, GRID, MUTED = '#1f2933', '#e3e7eb', '#7b8794'
MARK = ['o', 's', '^']
PAIR_COL = {'GATK_DV': '#999933', 'DV_PPR': '#44AA99', 'GATK_PPR': '#AA4499'}   # the pair colours of Figure 3
PIPE_COL = {'GATK': '#882255', 'DV': '#332288', 'PPR': '#117733'}               # the pipeline colours of Figures 1 and 4
PAIRS = [('GATK_PPR', 'GATK vs PPR'), ('DV_PPR', 'DV vs PPR'), ('GATK_DV', 'GATK vs DV')]
PIPES = ['GATK', 'DV', 'PPR']
TITLES = {'a': 'Depth on the two platforms', 'b': 'Concordance within a platform and depth',
          'c': 'Concordance across platforms', 'd': 'K-mer support and depth'}
STRATA = ['genome', 'BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']   # the order of Figure 1
STRATUM_LABEL = dict({s: s for s in STRATA}, genome='Genome')
LONG_DASHES = {0x2012, 0x2013, 0x2014, 0x2015, 0x2212}


def style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.0, 'axes.labelsize': 7.0, 'xtick.labelsize': 6.5,
                         'ytick.labelsize': 6.5, 'axes.edgecolor': '#9aa5b1', 'axes.labelcolor': INK, 'xtick.color': INK,
                         'ytick.color': INK, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
                         'axes.unicode_minus': False})


def read(name):
    with open(os.path.join(HERE, name), newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def load():
    """D: variable -> (individual, platform) -> value; COR: (family, x, y, platform, type, stratum) -> the correlation row;
    IND: the individuals (numbered 1-23)"""
    data = defaultdict(dict)
    for r in read('coverage_data.tsv'):
        data[r['variable']][(int(r['individual']), r['platform'])] = float(r['value'])
    rows = {(r['family'], r['x'], r['y'], r['platform'], r['type'], r['stratum']): r for r in read('coverage_correlations.tsv')}
    return data, rows, sorted({i for (i, _) in data['depth']})


D, COR, IND = load()


def depth(pl):
    return np.array([D['depth'][(i, pl)] for i in IND])


def cor(family, x, y, pl, typ='ALL', st='genome'):
    r = COR[(family, x, y, pl, typ, st)]
    return float(r['spearman_rho']), float(r['p_holm'])


def stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


def fit(ax, x, y, col):
    """Theil-Sen guide line (median of the pairwise slopes)"""
    i, j = np.triu_indices(len(x), 1)
    ok = x[j] != x[i]
    b = np.median((y[j][ok] - y[i][ok]) / (x[j][ok] - x[i][ok]))
    a = np.median(y - b * x)
    xx = np.linspace(x.min(), x.max(), 10)
    ax.plot(xx, a + b * xx, color=col, lw=1.0, alpha=0.85, zorder=2)


def finish(ax, xlabel, ylabel=None):
    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def series(ax, x, y, col, mk, label):
    ax.scatter(x, y, s=16, color=col, marker=mk, lw=0, alpha=0.95, zorder=3)
    fit(ax, x, y, col)
    return Line2D([], [], color=col, marker=mk, ls='-', lw=1.0, ms=4.0, label=label)


# ---------------------------------------------------------------- Figure 6
def draw_a(ax):
    """a: the mean depth of each individual, Aviti against Illumina, the line of equality, Spearman rho and p (one test)"""
    x, y = depth('illumina'), depth('aviti')
    lim = (0, max(x.max(), y.max()) * 1.08)
    ax.plot(lim, lim, color=MUTED, lw=1.0, ls='--')
    ax.scatter(x, y, s=16, color=INK, lw=0)
    rho, p = cor('depth', 'depth Illumina', 'depth Aviti', 'both', '')
    ax.text(0.97, 0.03, f'rho {rho:.2f}, p = {p:.2g}', transform=ax.transAxes, ha='right', va='bottom', fontsize=7, color=INK)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect('equal')
    finish(ax, 'Mean depth, Illumina (×)', 'Mean depth, Aviti (×)')


def draw_b(axes):
    """b: the genotype concordance of the three pipeline pairs within a platform against the depth of that platform; one axis per
    platform (Illumina, Aviti), the y axis shared"""
    for ax, (pl, name) in zip(axes, PLATS):
        x = depth(pl)
        hand = []
        for (pid, lab), mk in zip(PAIRS, MARK):
            y = np.array([D[f'concordance:{pid}:ALL:genome'][(i, pl)] for i in IND])
            rho, p = cor('within', f'depth {pl}', pid, pl)
            hand.append(series(ax, x, y, PAIR_COL[pid], mk, f'{lab}: rho {rho:.2f} {stars(p)}'))
        ax.legend(handles=hand, loc='lower right', fontsize=7, handlelength=1.6, borderaxespad=0.3)
        ax.set_title(name, fontsize=8)
        finish(ax, f'Mean depth, {name} (×)')
    axes[0].set_ylabel('Genotype concordance (%)')


def draw_c(ax):
    """c: the concordance of each pipeline between Illumina and Aviti against the lower of the two depths"""
    x = np.minimum(depth('illumina'), depth('aviti'))
    hand = []
    for s, mk in zip(PIPES, MARK):
        y = np.array([D[f'concordance:{s}_across:ALL:genome'][(i, 'both')] for i in IND])
        rho, p = cor('across', 'lower depth', f'{s} Illumina vs Aviti', 'both')
        hand.append(series(ax, x, y, PIPE_COL[s], mk, f'{s}: rho {rho:.2f} {stars(p)}'))
    ax.legend(handles=hand, loc='lower right', fontsize=7, handlelength=1.6, borderaxespad=0.3)
    finish(ax, 'Lower depth of the two platforms (×)', 'Concordance across platforms (%)')


def draw_d(axes):
    """d: the share of the records with k-mer support (Merfin) of the call set of each pipeline against depth; one axis per
    platform"""
    for ax, (pl, name) in zip(axes, PLATS):
        hand = []
        for s, mk in zip(PIPES, MARK):
            key = f'merfin:{s}_ALL'
            x = np.array([D['depth'][(i, pl)] for i in IND])
            y = np.array([D[key][(i, pl)] for i in IND])
            rho, p = cor('merfin', f'depth {pl}', f'{s}_ALL', pl)
            hand.append(series(ax, x, y, PIPE_COL[s], mk, f'{s}: rho {rho:.2f} {stars(p)}'))
        ax.legend(handles=hand, loc='lower right', fontsize=7, handlelength=1.4, borderaxespad=0.3)
        ax.set_title(name, fontsize=8)
        finish(ax, f'Mean depth, {name} (×)')
    axes[0].set_ylabel('Records with k-mer support (%)')


# ---------------------------------------------------------------- supplementary figure by region
def draw_regions(ax, rows, title):
    """Spearman rho by region: one point per series, filled when the Holm-adjusted p is below 0.05 (genome: Holm within the panel
    family of Figure 6; regions: Holm over the 10 regions), open otherwise; rows: (label, colour, marker, key of a stratum)"""
    offs = np.linspace(-0.22, 0.22, len(rows))
    hand = []
    for (label, col, mk, keyfn), dx in zip(rows, offs):
        for k, st in enumerate(STRATA):
            r = COR[keyfn(st)]
            rho, p = float(r['spearman_rho']), float(r['p_holm'])
            ax.scatter([k + dx], [rho], s=20, marker=mk, color=col if p < 0.05 else 'white', edgecolor=col, lw=1.0, zorder=3)
        hand.append(Line2D([], [], color=col, marker=mk, ls='', ms=4.2, label=label))
    ax.axhline(0, color=MUTED, lw=0.7)
    ax.axvline(0.5, color='#c3cbd3', lw=0.7)
    ax.set_xticks(range(len(STRATA)))
    ax.set_xticklabels([STRATUM_LABEL[s] for s in STRATA])
    ax.set_xlim(-0.6, len(STRATA) - 0.4)
    ax.set_ylim(-1.05, 1.05)
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax.set_ylabel('Spearman rho with depth')
    ax.grid(axis='y', color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(title, loc='left', fontsize=7.6, fontweight='bold')
    return hand


# ---------------------------------------------------------------- output
def no_long_dashes(fig):
    """stop if any text of the figure holds a long dash or the minus sign U+2212 (drawn first, so that the tick labels are set)"""
    fig.canvas.draw()
    bad = [t.get_text() for t in fig.findobj(Text) if any(ord(ch) in LONG_DASHES for ch in t.get_text())]
    if bad:
        sys.exit(f'long dash in the figure text: {bad[:3]}')


def save(fig, stem, tiff=False, tight=True):
    """figures/<stem>.pdf (vector) and .png at DPI dots per inch, with tiff also .tiff at DPI (LZW, RGB); after the text check;
    tight=False: the page keeps the size of the figure (drawn at its printed size)"""
    no_long_dashes(fig)
    os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True)
    kinds = [('pdf', {}), ('png', {'dpi': DPI})]
    if tiff:
        kinds.append(('tiff', {'dpi': DPI, 'pil_kwargs': {'compression': 'tiff_lzw'}}))
    for ext, kw in kinds:
        fig.savefig(os.path.join(HERE, 'figures', f'{stem}.{ext}'), bbox_inches='tight' if tight else None, pad_inches=0.04,
                    facecolor='white', **kw)
    if tiff:
        tif = os.path.join(HERE, 'figures', f'{stem}.tiff')
        Image.open(tif).convert('RGB').save(tif, compression='tiff_lzw', dpi=(DPI, DPI))   # flatten the alpha channel


def size_mm(fig):
    """the size of the saved figure in mm: the tight box plus the padding"""
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    return (bb.width + 0.08) * 25.4, (bb.height + 0.08) * 25.4
