#!/usr/bin/env python3
"""Figure 3: a genome-wide genotype concordance of the nine call-set pairs, one dot per individual; b same pipeline across platforms
minus same platform across pipelines, by region; c genotype concordance across platforms of each pipeline, by region (three styles).
Markers: in b circles for the series computed from GATK and PPR, squares for DV and PPR; in c circles GATK, squares DV, triangles PPR.
Inputs: concordance_data.tsv, concordance_tests.tsv (concordance_table.py). Usage: python3 figure3_panels.py [--c rows|columns|dots]
"""
import csv
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402
from matplotlib.transforms import offset_copy  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REGIONS = ['BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']   # the order of Figure 1
ROWS = ['genome'] + REGIONS
ROW_LABEL = dict({r: r for r in REGIONS}, genome='Genome')
TYPES = [('ALL', 'All variants'), ('SNV', 'SNVs'), ('INDEL', 'Indels')]
# Paul Tol's muted colours, as in Figure 1: the pipelines take the dark colours, the within-platform pairs three other colours
PIPE_COLOR = {'GATK': '#882255', 'DV': '#332288', 'PPR': '#117733'}
STEP_COLOR = {'GATK_DV': '#999933', 'DV_PPR': '#44AA99', 'GATK_PPR': '#AA4499'}
PIPE_MARKER = {'GATK': 'o', 'DV': 's', 'PPR': '^'}   # one shape per pipeline, as in Figure 6, so that colour is not the only cue
PANEL_TITLE = {'a': 'Genome-wide genotype concordance', 'b': 'Platform versus pipeline, by region',
               'c': 'Genotype concordance across platforms, by region'}
C_STYLES = ('rows', 'columns', 'dots')   # panel c: lines with the regions in rows (as b), lines along the x axis (as Figure 1), dots
STEP_LABEL = {'GATK_DV': 'GATK vs DV: the caller changes', 'DV_PPR': 'DV vs PPR: the reference changes',
              'GATK_PPR': 'GATK vs PPR: both change'}
POS, NEG, INK, GRID, NS = '#0072B2', '#D55E00', '#1f2933', '#e3e7eb', '#7b8794'
PIPES = ['GATK', 'DV', 'PPR']
# PPR against GATK and against DV: test id, the other pipeline (its colour), label
H1_SERIES = [('H1_crossplatform_vs_crossmethod', 'GATK', 'computed from GATK and PPR'),
             ('H1_DV_PPR_crossplatform_vs_crossmethod', 'DV', 'computed from DV and PPR')]
H2_MARKS = [('H2_PPR_vs_GATK_crossplatform', 'GATK'), ('H2_PPR_vs_DV_crossplatform', 'DV')]
STAR_X = (1.03, 1.29)      # the two columns of marks, in axes fraction (headers of 7 pt side by side)
# panel a: the nine pairs under three headers
PAIR_ROWS = [(None, 'Illumina'), ('GATK_DV_illumina', 'GATK vs DV'), ('DV_PPR_illumina', 'DV vs PPR'), ('GATK_PPR_illumina', 'GATK vs PPR'),
             (None, 'Aviti'), ('GATK_DV_aviti', 'GATK vs DV'), ('DV_PPR_aviti', 'DV vs PPR'), ('GATK_PPR_aviti', 'GATK vs PPR'),
             (None, 'Illumina vs Aviti'), ('GATK_illumina_aviti', 'GATK'), ('DV_illumina_aviti', 'DV'), ('PPR_illumina_aviti', 'PPR')]


def style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.titlesize': 8.5, 'axes.labelsize': 7.5,
                         'axes.edgecolor': '#9aa5b1', 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
                         'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
                         'axes.unicode_minus': False})


def load(folder=HERE):
    """D: (metric, pair, stratum, type) -> 23 values (fractions, individual 1-23); T: (comparison_id, metric, stratum, type) -> test row"""
    raw = {}
    with open(os.path.join(folder, 'concordance_data.tsv'), newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            for metric in ('genotype', 'allele'):
                raw.setdefault((metric, r['pair'], r['stratum'], r['type']), {})[int(r['individual'])] = float(r[f'{metric}_concordance'])
    D = {}
    for k, v in raw.items():
        if sorted(v) != list(range(1, 24)):
            raise SystemExit(f'{k}: {len(v)} individuals')
        D[k] = np.array([v[i] for i in range(1, 24)])
    with open(os.path.join(folder, 'concordance_tests.tsv'), newline='') as fh:
        T = {(r['comparison_id'], r['metric'], r['stratum'], r['type']): r for r in csv.DictReader(fh, delimiter='\t')}
    return D, T


def stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


def region_rows():
    """y positions of the genome row and the 10 regions (a gap below the genome row)"""
    ys = np.arange(len(ROWS))[::-1].astype(float)
    ys[1:] -= 0.6
    return dict(zip(ROWS, ys))


def row_axes(ax, ypos, labels=True):
    ys = [ypos[r] for r in ROWS]
    ax.set_yticks(ys)
    ax.set_yticklabels([ROW_LABEL[r] for r in ROWS])
    ax.tick_params(axis='y', length=0, labelleft=labels)
    ax.axhline(ypos['genome'] - 0.8, color='#c3cbd3', lw=1.0)
    ax.set_ylim(min(ys) - 0.7, max(ys) + 0.7)
    ax.grid(axis='x', color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def star_column(ax, T, cid, typ, ypos, metric='genotype', x=1.03, color=INK, dy=0.0, fs=6.5):
    for r in ROWS:
        p = float(T[(cid, metric, r, typ)]['p_holm'])
        ax.text(x, ypos[r] + dy, stars(p), transform=ax.get_yaxis_transform(), ha='left', va='center', fontsize=fs,
                fontweight='bold' if p < 0.05 else 'normal', color=color)


def star_headers(ax, others):
    """the pipeline compared with PPR, in its colour, over each column of marks"""
    for other, x in zip(others, STAR_X):
        ax.text(x, 1.0, other, transform=ax.transAxes, ha='left', va='bottom', fontsize=7, fontweight='bold', color=PIPE_COLOR[other])


def forest(ax, T, cid, typ, labels=True, metric='genotype'):
    """paired difference x - y per region (percentage points); filled = significant after Holm, hollow = not"""
    ypos = region_rows()
    for r in ROWS:
        t = T[(cid, metric, r, typ)]
        md, p = float(t['mean_diff_pp']), float(t['p_holm'])
        c = POS if md >= 0 else NEG
        ax.hlines(ypos[r], 0, md, color=c, lw=1.4)
        ax.plot(md, ypos[r], 'o', ms=4.6, color=c, mfc=c if p < 0.05 else 'white', mew=1.1, zorder=3)
    ax.axvline(0, color=NS, lw=0.7, ls=(0, (3, 2)))
    row_axes(ax, ypos, labels)
    star_column(ax, T, cid, typ, ypos, metric)


def forest_series(ax, T, series, typ, labels=True, metric='genotype'):
    """several paired differences per region, one colour per series (the other pipeline compared with PPR); filled = significant
    after Holm, hollow = not; one column of marks per series, in its colour"""
    ypos = region_rows()
    offs = np.linspace(0.18, -0.18, len(series))
    for (cid, other, _), dy, sx in zip(series, offs, STAR_X):
        c = PIPE_COLOR[other]
        for r in ROWS:
            t = T[(cid, metric, r, typ)]
            md, p = float(t['mean_diff_pp']), float(t['p_holm'])
            ax.hlines(ypos[r] + dy, 0, md, color=c, lw=1.3)
            ax.plot(md, ypos[r] + dy, PIPE_MARKER[other], ms=4.0, color=c, mfc=c if p < 0.05 else 'white', mew=1.0, zorder=3)
        star_column(ax, T, cid, typ, ypos, metric, x=sx, color=c, dy=dy, fs=7)
    ax.axvline(0, color=NS, lw=1.0, ls=(0, (3, 2)))
    row_axes(ax, ypos, labels)
    star_headers(ax, [other for _, other, _ in series])


def draw_pairs(ax, D, typ, labels=True):
    """panel a: genome-wide genotype concordance (%) of the nine pairs, one dot per individual, bar = mean"""
    y, ys, ticks = 0.0, [], []
    for k, (pair, lab) in enumerate(PAIR_ROWS):
        if pair is None and k:
            y -= 0.35
        ys.append(y)
        ticks.append(lab)
        y -= 1
    jit = (np.arange(23) % 5 - 2) * 0.07
    for yy, (pair, lab) in zip(ys, PAIR_ROWS):
        if pair is None:
            continue
        v = 100 * D[('genotype', pair, 'genome', typ)]
        if pair.endswith('illumina_aviti'):
            c, mk, mfc = PIPE_COLOR[pair.split('_')[0]], 'D', None
        else:
            plat = pair.rsplit('_', 1)[1]
            c = STEP_COLOR['_'.join(pair.split('_')[:2])]
            mk, mfc = ('o', None) if plat == 'illumina' else ('s', 'white')
        ax.plot(v, yy + jit, ls='', marker=mk, ms=3.4, color=c, mfc=mfc or c, mew=1.0, alpha=0.9)
        ax.plot([v.mean()] * 2, [yy - 0.36, yy + 0.36], color=INK, lw=1.4, zorder=4)
    ax.set_yticks(ys)
    ax.set_yticklabels(ticks)
    for t, (pair, _) in zip(ax.get_yticklabels(), PAIR_ROWS):
        if pair is None:
            t.set_fontweight('bold')
    ax.tick_params(axis='y', length=0, labelleft=labels)
    ax.set_ylim(min(ys) - 0.7, max(ys) + 0.6)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
    ax.grid(axis='x', color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def draw_same_pipeline_vs_platform(ax, T, typ, labels=True):
    """panel b: same pipeline across platforms minus same platform across pipelines, computed from GATK and PPR and from DV and PPR"""
    forest_series(ax, T, H1_SERIES, typ, labels)


def cross_means(D, typ):
    """genotype concordance across platforms (%) of each pipeline per row, mean of the 23 individuals"""
    return {s: {r: 100 * D[('genotype', f'{s}_illumina_aviti', r, typ)].mean() for r in ROWS} for s in PIPES}


def h2_marks_rows(ax, T, typ, ypos):
    for (cid, other), sx in zip(H2_MARKS, STAR_X):
        star_column(ax, T, cid, typ, ypos, x=sx, color=PIPE_COLOR[other])
    star_headers(ax, [other for _, other in H2_MARKS])


def draw_across_platforms_dots(ax, D, T, typ, labels=True):
    """panel c (dots): genotype concordance across platforms (%) of each pipeline, mean of the 23 individuals; marks: PPR vs GATK
    (GATK colour) and PPR vs DV (DV colour)"""
    ypos = region_rows()
    off = {'GATK': 0.2, 'DV': 0.0, 'PPR': -0.2}
    m = cross_means(D, typ)
    for r in ROWS:
        for s in PIPES:
            ax.plot(m[s][r], ypos[r] + off[s], PIPE_MARKER[s], ms=3.8, color=PIPE_COLOR[s], zorder=3)
    row_axes(ax, ypos, labels)
    h2_marks_rows(ax, T, typ, ypos)


def draw_across_platforms_rows(ax, D, T, typ, labels=True):
    """panel c (lines, the regions in rows as in b): the 10 regions of each pipeline joined by a line, the genome apart; marks: PPR vs
    GATK and PPR vs DV in columns on the right"""
    ypos = region_rows()
    off = {'GATK': 0.15, 'DV': 0.0, 'PPR': -0.15}
    m = cross_means(D, typ)
    for s in PIPES:
        ax.plot([m[s][r] for r in REGIONS], [ypos[r] + off[s] for r in REGIONS], '-', color=PIPE_COLOR[s], lw=1.0, alpha=0.9, zorder=2)
        ax.plot([m[s][r] for r in ROWS], [ypos[r] + off[s] for r in ROWS], ls='', marker=PIPE_MARKER[s], ms=3.4,
                color=PIPE_COLOR[s], zorder=3)
    row_axes(ax, ypos, labels)
    h2_marks_rows(ax, T, typ, ypos)


def draw_across_platforms_columns(ax, D, T, typ, labels=True):
    """panel c (lines along the x axis, as Figure 1): the 10 regions of each pipeline joined by a line, the genome apart on the left;
    marks: PPR vs GATK (upper row) and PPR vs DV (lower row) above the plot"""
    xs = dict({'genome': 0.0}, **{r: 1.3 + i for i, r in enumerate(REGIONS)})
    m = cross_means(D, typ)
    for s in PIPES:
        ax.plot([xs[r] for r in REGIONS], [m[s][r] for r in REGIONS], '-', color=PIPE_COLOR[s], lw=1.0, alpha=0.9, zorder=2)
        ax.plot([xs[r] for r in ROWS], [m[s][r] for r in ROWS], ls='', marker=PIPE_MARKER[s], ms=4.0, color=PIPE_COLOR[s], zorder=3)
    ax.axvline(0.65, color='#c3cbd3', lw=1.0)
    ax.set_xticks([xs[r] for r in ROWS])
    ax.set_xticklabels([ROW_LABEL[r] for r in ROWS], rotation=90)
    ax.tick_params(axis='x', length=0)
    ax.tick_params(axis='y', labelleft=labels)
    ax.set_xlim(-0.6, xs[REGIONS[-1]] + 0.6)
    ax.grid(axis='y', color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for k, (cid, other) in enumerate(H2_MARKS):
        dy = 3 + 10 * (len(H2_MARKS) - 1 - k)      # points above the axes: rows of 7-pt marks 10 pt apart, whatever the axes height
        trans = offset_copy(ax.get_xaxis_transform(), ax.figure, y=dy, units='points')
        for r in ROWS:
            p = float(T[(cid, 'genotype', r, typ)]['p_holm'])
            ax.text(xs[r], 1.0, stars(p), transform=trans, ha='center', va='bottom', fontsize=7, color=PIPE_COLOR[other],
                    fontweight='bold' if p < 0.05 else 'normal')
        if labels:
            ax.text(-0.04, 1.0, f'vs {other}', transform=offset_copy(ax.transAxes, ax.figure, y=dy, units='points'), ha='right',
                    va='bottom', fontsize=7, fontweight='bold', color=PIPE_COLOR[other])


def draw_across_platforms(ax, D, T, typ, style, labels=True):
    {'rows': draw_across_platforms_rows, 'columns': draw_across_platforms_columns,
     'dots': draw_across_platforms_dots}[style](ax, D, T, typ, labels)


def c_style(argv=None):
    """--c rows|columns|dots from the command line; columns (lines as in Figure 1) by default"""
    argv = sys.argv if argv is None else argv
    style = argv[argv.index('--c') + 1] if '--c' in argv else 'columns'
    if style not in C_STYLES:
        raise SystemExit(f'--c must be one of {", ".join(C_STYLES)}')
    return style


def legend_pairs():
    h = [Line2D([], [], ls='', marker='o', ms=4.5, color=STEP_COLOR[s], label=STEP_LABEL[s]) for s in STEP_COLOR]
    h += [Line2D([], [], ls='', marker='D', ms=4.5, color=PIPE_COLOR[s], label=f'{s}: Illumina vs Aviti') for s in PIPES]
    h += [Line2D([], [], ls='', marker='o', ms=4.5, color='#555555', label='Illumina pair'),
          Line2D([], [], ls='', marker='s', ms=4.5, color='#555555', mfc='white', label='Aviti pair'),
          Line2D([], [], color=INK, lw=1.4, label='mean of the 23 individuals')]
    return h


def legend_forest(pos_label, neg_label):
    """colour keys only: the significance keys (open points, the marks) are in the figure legend of the manuscript"""
    return [Line2D([], [], ls='', marker='o', ms=5, color=POS, label=pos_label),
            Line2D([], [], ls='', marker='o', ms=5, color=NEG, label=neg_label)]


def legend_series():
    """panel b: the colour of each series (the significance keys are in the figure legend of the manuscript)"""
    return [Line2D([], [], ls='-', lw=1.3, marker=PIPE_MARKER[other], ms=4.5, color=PIPE_COLOR[other], label=lab)
            for _, other, lab in H1_SERIES]


def legend_pipelines(lines=False):
    return [Line2D([], [], ls='-' if lines else '', lw=1.0, marker=PIPE_MARKER[s], ms=4.5, color=PIPE_COLOR[s],
                   label=f'{s}: Illumina vs Aviti')
            for s in PIPES]


def save(fig, stem):
    os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True)
    for ext, kw in (('pdf', {}), ('png', {'dpi': 600})):
        fig.savefig(os.path.join(HERE, 'figures', f'{stem}.{ext}'), bbox_inches='tight', facecolor='white', **kw)
    plt.close(fig)
    print('written', stem)


def panel_title(fig, tag):
    fig.suptitle(PANEL_TITLE[tag], x=0.02, ha='left', fontweight='bold', fontsize=8.5, y=1.04)


if __name__ == '__main__':
    style()
    D, T = load()
    cs = c_style()
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.6), sharey=True)
    for ax, (typ, tl) in zip(axes, TYPES):
        draw_pairs(ax, D, typ, ax is axes[0])
        ax.set_title(tl, fontweight='bold')
    panel_title(fig, 'a')
    fig.supxlabel('Genotype concordance in the genome (%)', fontsize=7.5, y=-0.02)
    fig.legend(handles=legend_pairs(), loc='upper center', bbox_to_anchor=(0.5, -0.06), ncol=3, fontsize=6.3)
    save(fig, 'Figure3a_pairs')
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.9), sharey=True)
    for ax, (typ, tl) in zip(axes, TYPES):
        draw_same_pipeline_vs_platform(ax, T, typ, ax is axes[0])
        ax.set_title(tl, fontweight='bold', pad=10)
    panel_title(fig, 'b')
    fig.subplots_adjust(wspace=0.45, right=0.9)
    fig.supxlabel('Same pipeline across platforms minus same platform across pipelines (percentage points)', fontsize=7.5, y=-0.02)
    fig.legend(handles=legend_series(), loc='upper center', bbox_to_anchor=(0.5, -0.06), ncol=2, fontsize=6.3)
    save(fig, 'Figure3b_same_pipeline_vs_same_platform')
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 3.0 if cs == 'columns' else 2.7), sharey=True)
    for ax, (typ, tl) in zip(axes, TYPES):
        draw_across_platforms(ax, D, T, typ, cs, ax is axes[0])
        ax.set_title(tl, fontweight='bold', pad=17 if cs == 'columns' else 10)
    panel_title(fig, 'c')
    if cs == 'columns':
        fig.subplots_adjust(wspace=0.12)
        axes[0].set_ylabel('Genotype concordance across platforms (%)')
    else:
        fig.subplots_adjust(wspace=0.45, right=0.9)
        fig.supxlabel('Genotype concordance across platforms (%)', fontsize=7.5, y=-0.02)
    fig.legend(handles=legend_pipelines(cs != 'dots'), loc='upper center', bbox_to_anchor=(0.5, -0.06 if cs != 'columns' else -0.12),
               ncol=3, fontsize=6.3)
    save(fig, 'Figure3c_across_platforms')
    print('FIGURE3_PANELS_DONE')
