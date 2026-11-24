#!/usr/bin/env python3
"""Figure 1a (SNVs) and 1b (indels): precision, recall and F1 of GATK, BWA-MEM + DeepVariant and PPR on HG002, Illumina and Aviti, per
region. PASS calls as solid lines with filled points; GATK before VQSR (all calls) as dashed lines with open points.
Inputs: benchmark_data.tsv. Outputs: figures/Figure1a_SNV.{pdf,png}, figures/Figure1b_INDEL.{pdf,png} (PNG 600 dpi).
Usage: python3 figure1_panels.py   (figure1_combine.py imports load() and draw_column() from this file)
"""
import csv
import math
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator, PercentFormatter  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REGIONS = ['BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']
PIPES = ['GATK', 'BWA-MEM + DeepVariant', 'PPR']
PLATFORMS = ['Illumina', 'Aviti']
# Paul Tol's muted colour-blind-safe scheme: a dark (Illumina) and a light (Aviti) colour per pipeline; circles (Illumina) and
# squares (Aviti) as in Supplementary Figures S1 and S2
COLOR = {('GATK', 'Illumina'): '#882255', ('GATK', 'Aviti'): '#CC6677',
         ('BWA-MEM + DeepVariant', 'Illumina'): '#332288', ('BWA-MEM + DeepVariant', 'Aviti'): '#88CCEE',
         ('PPR', 'Illumina'): '#117733', ('PPR', 'Aviti'): '#44AA99'}
MARKER = {'Illumina': 'o', 'Aviti': 's'}
MSIZE = {'Illumina': 4.5, 'Aviti': 4.0}   # a square looks larger than a circle of the same size
METRICS = [('precision', 'Precision'), ('recall', 'Recall'), ('f1', 'F1')]
TYPE_LABEL = {'SNP': 'SNVs', 'INDEL': 'Indels'}
INK, GRID = '#1f2933', '#e3e7eb'


def style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 8,
                         'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7, 'axes.linewidth': 1.0,
                         'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
                         'axes.edgecolor': '#9aa5b1', 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
                         'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
                         'axes.unicode_minus': False})


def load(path=os.path.join(HERE, 'benchmark_data.tsv')):
    """(platform, pipeline, region, type, filter) -> {precision, recall, f1} as percentages"""
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    if len(rows) != 160 or sorted({r['region'] for r in rows}) != sorted(REGIONS):
        raise SystemExit(f'unexpected table {path}: {len(rows)} rows')
    return {(r['platform'], r['pipeline'], r['region'], r['type'], r['filter']): {m: 100 * float(r[m]) for m, _ in METRICS}
            for r in rows}


def series(data, vtype, metric):
    """the lines of one panel: (platform, pipeline, filter, values in region order)"""
    out = []
    for plat in PLATFORMS:
        for pipe in PIPES:
            out.append((plat, pipe, 'PASS', [data[(plat, pipe, r, vtype, 'PASS')][metric] for r in REGIONS]))
        out.append((plat, 'GATK', 'ALL', [data[(plat, 'GATK', r, vtype, 'ALL')][metric] for r in REGIONS]))
    return out


def y_limits(values):
    """lower end: the lowest value minus 1 point, down to a multiple of 5 %; upper end 100 % plus a margin"""
    lo = 5 * math.floor((min(values) - 1) / 5)
    return lo, 100.8


def draw_column(axes, data, vtype, show_ylabel=True, show_xlabels=True):
    """one variant type in three stacked axes (precision, recall, F1)"""
    x = np.arange(len(REGIONS))
    for ax, (metric, mlab) in zip(axes, METRICS):
        vals = []
        for plat, pipe, filt, y in series(data, vtype, metric):
            vals += y
            c = COLOR[(pipe, plat)]
            if filt == 'PASS':
                ax.plot(x, y, color=c, lw=1.3, ls='-', marker=MARKER[plat], ms=MSIZE[plat], zorder=3)
            else:
                ax.plot(x, y, color=c, lw=1.1, ls=(0, (3, 2)), marker=MARKER[plat], ms=MSIZE[plat] + 0.2, mfc='white', mew=1.0,
                        zorder=2)
        lo, hi = y_limits(vals)
        ax.set_ylim(lo, hi)
        span = 100 - lo
        ax.yaxis.set_major_locator(MultipleLocator(5 if span > 15 else 2))
        ax.yaxis.set_major_formatter(PercentFormatter(100, decimals=0))
        ax.grid(color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.set_xlim(-0.5, len(REGIONS) - 0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(REGIONS if (show_xlabels and ax is axes[-1]) else [], rotation=45, ha='right', rotation_mode='anchor')
        if show_ylabel:
            ax.set_ylabel(mlab, fontweight='bold')
    axes[0].set_title(TYPE_LABEL[vtype], fontweight='bold')


def legend_handles():
    h = [Line2D([], [], color=COLOR[(pipe, plat)], lw=1.4, marker=MARKER[plat], ms=MSIZE[plat] + 0.5, label=f'{pipe}, {plat}')
         for pipe in PIPES for plat in PLATFORMS]
    h.append(Line2D([], [], color='#444444', lw=1.3, marker='o', ms=MSIZE['Illumina'], label='PASS calls'))
    h.append(Line2D([], [], color='#444444', lw=1.1, ls=(0, (3, 2)), marker='o', ms=MSIZE['Illumina'], mfc='white', mew=1.0,
                    label='GATK before VQSR (all calls)'))
    return h


def save(fig, stem):
    os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True)
    for ext, kw in (('pdf', {}), ('png', {'dpi': 600})):
        fig.savefig(os.path.join(HERE, 'figures', f'{stem}.{ext}'), bbox_inches='tight', facecolor='white', **kw)
    print('written', stem)


if __name__ == '__main__':
    style()
    d = load()
    for vtype, stem in (('SNP', 'Figure1a_SNV'), ('INDEL', 'Figure1b_INDEL')):
        fig, axes = plt.subplots(3, 1, figsize=(3.6, 6.6), sharex=True)
        draw_column(axes, d, vtype)
        fig.legend(handles=legend_handles(), loc='lower center', ncol=2, bbox_to_anchor=(0.5, -0.12), fontsize=6.5)
        fig.subplots_adjust(hspace=0.12)
        save(fig, stem)
        plt.close(fig)
    print('FIGURE1_PANELS_DONE')
