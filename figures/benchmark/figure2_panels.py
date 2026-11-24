#!/usr/bin/env python3
"""The panels of Figure 2: a, b precision-recall curves (SNVs, indels); c, d errors in the benchmark regions; e reads mapped and reads with
MAPQ 60 of the five set-ups in the BAM files that DeepVariant read, on chr1-22, X and Y: the giraffe alignments after vg surject and the
BWA-MEM alignments on the same chromosomes (all metrics in Supplementary Figure S7). Every point is a filled circle in the colour of its
set-up; Illumina solid lines, Aviti dashed lines; in c and d the error kinds differ by the colour of the line; in e a stem in the line
style of the platform leads to each point.
Inputs: graph_arms_data.tsv, graph_arms_roc.tsv, graph_arms_mapping.tsv. Outputs: figures/Figure2a_PR_SNV ... Figure2e_alignment (.pdf, .png).
Usage: python3 figure2_panels.py   (figure2_combine.py and supp_graph_arms_alignment.py import this file)
"""
import csv
import math
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MaxNLocator, MultipleLocator  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import graph_setups as gs  # noqa: E402

TYPE_LABEL = {'SNP': 'SNVs', 'INDEL': 'Indels'}
INK, GRID = '#1f2933', '#e3e7eb'
LS = {'Illumina': '-', 'Aviti': (0, (3, 2))}
# c, d: one line colour per error kind, darker than the set-up colours (Okabe-Ito) of the points; the total line is the thickest
ERR_COLOR = {'total': '#222222', 'FN': '#B2182B', 'FP': '#5E3C99'}
ERR_LW = {'total': 1.8, 'FN': 1.2, 'FP': 1.2}
ERR_LABEL = {'total': 'Total errors (false negatives + false positives)', 'FN': 'False negatives', 'FP': 'False positives'}
METRICS = [('gam_aligned', 'GAM aligned'), ('gam_perfect', 'GAM perfect'), ('gam_gapless', 'GAM gapless'),
           ('gam_properly_paired', 'GAM properly paired'), ('gam_mapq60', 'GAM MAPQ 60'), ('bam_mapped', 'BAM mapped'),
           ('bam_properly_paired', 'BAM properly paired'), ('bam_mapq60', 'BAM MAPQ 60')]
# BWA-MEM has only BAM metrics (no graph alignment): bam_* on chr1-22, X, Y (like the giraffe BAMs after surjection), bam_all_* on all
# contigs
BWA_METRICS = ['bam_mapped', 'bam_properly_paired', 'bam_mapq60', 'bam_all_mapped', 'bam_all_mapq60']
# Figure 2e: the BAM files that DeepVariant read, chr1-22, X, Y, for all five set-ups
FIG2_METRICS = ['bam_mapped', 'bam_mapq60']
FIG2_TITLES = {'bam_mapped': 'Mapped', 'bam_mapq60': 'MAPQ 60'}
# Supplementary Figure S7: BWA-MEM in the GAM aligned and GAM MAPQ 60 panels, from its BAM on all contigs
BWA_IN_GAM = {'gam_aligned': 'bam_all_mapped', 'gam_mapq60': 'bam_all_mapq60'}


def style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 8,
                         'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7, 'axes.linewidth': 1.0,
                         'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
                         'axes.edgecolor': '#9aa5b1', 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
                         'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
                         'axes.unicode_minus': False})


def load():
    """{'bb': (platform, set-up, type) -> PASS values in the benchmark regions, 'roc': (platform, set-up, type) -> (recall, precision)
    sorted by recall, 'align': (platform, set-up, metric) -> % of reads}"""
    def read(name):
        with open(os.path.join(HERE, name), newline='') as fh:
            return list(csv.DictReader(fh, delimiter='\t'))

    data, roc, align = read('graph_arms_data.tsv'), read('graph_arms_roc.tsv'), read('graph_arms_mapping.tsv')
    n_align = len(gs.PLATFORMS) * (len(gs.MAPPING_SETUPS) * len(METRICS) + len(BWA_METRICS))
    if len(data) != 200 or len(align) != n_align:
        raise SystemExit(f'unexpected tables: {len(data)} rows in graph_arms_data.tsv, {len(align)} in graph_arms_mapping.tsv')
    bb = {(r['platform'], r['setup'], r['type']): {k: float(r[k]) for k in ('precision', 'recall', 'FN', 'FP')}
          for r in data if r['region'] == 'BB'}
    curves = {}
    for r in roc:
        curves.setdefault((r['platform'], r['setup'], r['type']), []).append((float(r['recall']), float(r['precision'])))
    if len(bb) != 20 or len(curves) != 20:
        raise SystemExit(f'unexpected tables: {len(bb)} PASS points, {len(curves)} curves (20 each expected)')
    curves = {k: tuple(map(list, zip(*sorted(v)))) for k, v in curves.items()}
    align_d = {(r['platform'], r['setup'], r['metric']): float(r['percent']) for r in align}
    for p in gs.PLATFORMS:
        for s in gs.SETUPS:
            for metric in FIG2_METRICS:
                if (p, s, metric) not in align_d:
                    raise SystemExit(f'graph_arms_mapping.tsv: no {metric} for {p} {s}')
    return {'bb': bb, 'roc': curves, 'align': align_d}


def marker_kw(s, ms):
    """every point a filled circle in the colour of its set-up, on both platforms (the platform is in the line style)"""
    c = gs.COLOR[s]
    return dict(marker='o', ms=ms, color=c, mfc=c, mec='white', mew=0.4)


def draw_pr(ax, d, vtype, show_ylabel=True):
    """precision against recall as the QUAL threshold rises (hap.py ROC, all records), and the point of the PASS calls"""
    for plat in gs.PLATFORMS:
        for s in gs.SETUPS:
            r, p = d['roc'][(plat, s, vtype)]
            ax.plot(r, p, color=gs.COLOR[s], lw=1.1, ls=LS[plat], zorder=2)
            b = d['bb'][(plat, s, vtype)]
            ax.plot(b['recall'], b['precision'], ls='', zorder=4, **marker_kw(s, 5.5))
    rec = [v['recall'] for k, v in d['bb'].items() if k[2] == vtype]
    prec = [v['precision'] for k, v in d['bb'].items() if k[2] == vtype]
    # the axes start below the lowest PASS point and end just after the highest, so the ends of the curves fill the panel
    ax.set_xlim(math.floor((min(rec) - 0.025) / 0.01) * 0.01, min(1.0, math.ceil((max(rec) + 0.002) / 0.005) * 0.005))
    ax.set_ylim(math.floor((min(prec) - 0.002) / 0.005) * 0.005, 1.0004)
    ax.xaxis.set_major_locator(MultipleLocator(0.01 if vtype == 'SNP' else 0.02))
    ax.yaxis.set_major_locator(MultipleLocator(0.002 if vtype == 'SNP' else 0.005))
    ax.set_xlabel('Recall')
    if show_ylabel:
        ax.set_ylabel('Precision')
    ax.set_title(TYPE_LABEL[vtype], fontweight='bold')
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def draw_errors(ax, d, vtype, show_ylabel=True):
    """false negatives, false positives and their sum (thousands) in the benchmark regions, PASS calls, set-ups in decomposition order;
    the error kind in the colour of the line, the set-up in the colour of the point"""
    x = np.arange(len(gs.SETUPS))
    top = 0
    for plat in gs.PLATFORMS:
        fn = np.array([d['bb'][(plat, s, vtype)]['FN'] for s in gs.SETUPS]) / 1000
        fp = np.array([d['bb'][(plat, s, vtype)]['FP'] for s in gs.SETUPS]) / 1000
        tot = fn + fp
        top = max(top, tot.max())
        for kind, vals in (('total', tot), ('FN', fn), ('FP', fp)):
            ax.plot(x, vals, color=ERR_COLOR[kind], lw=ERR_LW[kind], ls=LS[plat], zorder=3)
            for k, s in enumerate(gs.SETUPS):
                ax.plot(k, vals[k], ls='', zorder=4, **marker_kw(s, 5.6 if kind == 'total' else 4.8))
    ax.set_ylim(0, top * 1.08)
    ax.set_xlim(-0.35, len(x) - 0.65)
    ax.set_xticks(x)
    ax.set_xticklabels([gs.SHORT[s] for s in gs.SETUPS], rotation=30, ha='right', rotation_mode='anchor')
    if show_ylabel:
        ax.set_ylabel('Errors in the benchmark\nregions (thousands)')
    ax.set_title(TYPE_LABEL[vtype], fontweight='bold')
    ax.grid(axis='y', color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def draw_alignment(axes, d, metrics=FIG2_METRICS, ncols=4, setups=None, titles=None, bwa_in_gam=False):
    """one small panel per alignment metric (% of reads): the set-ups in rows (default all five), Illumina above Aviti in each row.
    titles: metric -> panel title (default the names of METRICS). A set-up without the metric (BWA-MEM in the GAM panels) gets 'no graph
    alignment' in its row; with bwa_in_gam, the GAM aligned and GAM MAPQ 60 panels show BWA-MEM from its BAM on all contigs instead,
    marked 'BAM, all contigs'"""
    setups = setups or gs.SETUPS
    titles = titles or dict(METRICS)
    ys = np.arange(len(setups))[::-1]
    off = {'Illumina': 0.17, 'Aviti': -0.17}

    def key(p, s, metric):
        if (p, s, metric) in d['align']:
            return (p, s, metric), False
        if bwa_in_gam and s == 'bwa_grch38' and metric in BWA_IN_GAM:
            return (p, s, BWA_IN_GAM[metric]), True
        return None, False

    for i, (ax, metric) in enumerate(zip(axes, metrics)):
        vals = [d['align'][key(p, s, metric)[0]] for p in gs.PLATFORMS for s in setups if key(p, s, metric)[0]]
        lo, hi = min(vals), max(vals)
        left, right = lo - (hi - lo) * 0.2 - 0.02, hi + (hi - lo) * 0.12 + 0.02
        for y, s in zip(ys, setups):
            k, substitute = key(gs.PLATFORMS[0], s, metric)
            if k is None:
                ax.text(left + (right - left) * 0.03, y, 'no graph alignment', fontsize=7, color='#7b8794', va='center')
                continue
            for plat in gs.PLATFORMS:
                v = d['align'][key(plat, s, metric)[0]]
                ax.hlines(y + off[plat], left, v, color=gs.COLOR[s], lw=1.0, alpha=0.5, ls=LS[plat])
                ax.plot(v, y + off[plat], ls='', zorder=3, **marker_kw(s, 5))
            if substitute:     # 6 pt and under the points, so that it does not cover the BWA-MEM point
                ax.text(left + (right - left) * 0.03, y, 'BAM, all contigs', fontsize=6, color='#52606d', va='center', zorder=2,
                        bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
        ax.set_xlim(left, right)
        ax.set_ylim(-0.55, len(ys) - 0.45)
        ax.set_yticks(ys)
        ax.set_yticklabels([gs.SHORT[s] for s in setups] if i % ncols == 0 else [])
        ax.tick_params(axis='y', length=0)
        ax.tick_params(axis='x', labelsize=7)
        ax.xaxis.set_major_locator(MaxNLocator(3))
        ax.set_title(f'{titles[metric]} (%)', fontweight='bold', fontsize=8)
        ax.grid(axis='x', color=GRID, lw=0.6)
        ax.set_axisbelow(True)


def setup_handles(setups=None):
    return [Line2D([], [], ls='', label=gs.LABEL[s], **marker_kw(s, 5.5)) for s in (setups or gs.SETUPS)]


def key_handles(kinds=('platform', 'errors')):
    """the keys that are not set-ups: the line style of each platform and the line colour of each error kind of c and d"""
    h = []
    if 'platform' in kinds:
        h += [Line2D([], [], color='#555555', lw=1.2, ls=LS['Illumina'], label='Illumina: solid lines'),
              Line2D([], [], color='#555555', lw=1.2, ls=LS['Aviti'], label='Aviti: dashed lines')]
    if 'platform' in kinds and 'errors' in kinds:
        h.append(Line2D([], [], ls='', label=' '))     # empty slot, so that the error kinds fill the second column
    if 'errors' in kinds:
        h += [Line2D([], [], color=ERR_COLOR[k], lw=ERR_LW[k], label=ERR_LABEL[k]) for k in ('total', 'FN', 'FP')]
    return h


def save(fig, stem):
    os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True)
    for ext, kw in (('pdf', {}), ('png', {'dpi': 600})):
        fig.savefig(os.path.join(HERE, 'figures', f'{stem}.{ext}'), bbox_inches='tight', facecolor='white', **kw)
    print('written', stem)


if __name__ == '__main__':
    style()
    d = load()
    for vtype, stem in (('SNP', 'Figure2a_PR_SNV'), ('INDEL', 'Figure2b_PR_INDEL')):
        fig, ax = plt.subplots(figsize=(3.5, 3.0))
        draw_pr(ax, d, vtype)
        fig.legend(handles=setup_handles() + key_handles(('platform',)), loc='upper center', bbox_to_anchor=(0.5, 0.0),
                   fontsize=5.5, title=gs.LEGEND_TITLE, title_fontsize=5.5, alignment='left')
        save(fig, stem)
        plt.close(fig)
    for vtype, stem in (('SNP', 'Figure2c_errors_SNV'), ('INDEL', 'Figure2d_errors_INDEL')):
        fig, ax = plt.subplots(figsize=(3.5, 3.0))
        draw_errors(ax, d, vtype)
        fig.legend(handles=setup_handles() + key_handles(('platform', 'errors')), loc='upper center', bbox_to_anchor=(0.5, -0.02),
                   fontsize=5.5, title=gs.LEGEND_TITLE, title_fontsize=5.5, alignment='left')
        save(fig, stem)
        plt.close(fig)
    fig, axes = plt.subplots(1, len(FIG2_METRICS), figsize=(7.1, 1.75))
    draw_alignment(axes, d, titles=FIG2_TITLES)
    fig.subplots_adjust(wspace=0.12)
    fig.legend(handles=setup_handles() + key_handles(('platform',)), loc='upper center', bbox_to_anchor=(0.5, -0.12),
               ncol=2, fontsize=6, alignment='left')
    save(fig, 'Figure2e_alignment')
    plt.close(fig)
    print('FIGURE2_PANELS_DONE')
