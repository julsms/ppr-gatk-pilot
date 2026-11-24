#!/usr/bin/env python3
"""genomescope_panels.py - style, data and drawing functions of the two GenomeScope supplementary figures (imported), with the tables of
genomescope_stats.py from this folder: paired lines per individual with the platform medians and the Holm-corrected tests, and the
heterozygosity against the genome length estimate. Platforms in blue (Illumina) and orange (Aviti).
"""
import csv
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DPI = int(os.environ.get('DPI', '600'))          # resolution of the PNG files
PLAT = {'illumina': ('Illumina', '#4c72b0'), 'aviti': ('Aviti', '#dd8452')}
EST = [('heterozygosity', 'Heterozygosity (%)'), ('error_rate', 'Read error rate (%)'), ('model_fit', 'Model fit (%)'),
       ('genome_length_gb', 'Genome length (Gb)')]
INK, GRID = '#1f2933', '#e3e7eb'


def style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.labelsize': 7.5, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                         'axes.edgecolor': '#9aa5b1', 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
                         'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'axes.unicode_minus': False})


def read(name, folder=HERE):
    with open(os.path.join(folder, name), newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def load(folder=HERE):
    """V: (platform, individual) -> the estimates of that individual on that platform; N: the number of individuals; ST: estimate -> its
    row of genomescope_stats.tsv (medians, W, p, Holm p, rank-biserial r); COR: the rows of genomescope_correlations.tsv"""
    ind = read('genomescope_per_individual.tsv', folder)
    V = {(r['platform'], int(r['individual'])): r for r in ind}
    N = max(int(r['individual']) for r in ind)
    if len(ind) != len(V) or sorted(V) != sorted((p, i) for p in PLAT for i in range(1, N + 1)):
        raise SystemExit('genomescope_per_individual.tsv: not one row per individual (1 to N) and platform')
    ST = {r['estimate']: r for r in read('genomescope_stats.tsv', folder)}
    for e, _ in EST:
        if int(ST[e]['n_pairs']) != N:
            raise SystemExit(f'genomescope_stats.tsv: {e} tested on {ST[e]["n_pairs"]} pairs, {N} individuals drawn')
    COR = read('genomescope_correlations.tsv', folder)
    return V, N, ST, COR


def stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


def save(fig, stem, folder=HERE):
    """figures/<stem>.pdf (vector) and figures/<stem>.png (DPI)"""
    os.makedirs(os.path.join(folder, 'figures'), exist_ok=True)
    for ext, kw in (('pdf', {}), ('png', {'dpi': DPI})):
        fig.savefig(os.path.join(folder, 'figures', f'{stem}.{ext}'), bbox_inches='tight', facecolor='white', **kw)
    plt.close(fig)
    print(f'written figures/{stem}.pdf and .png ({DPI} dpi)')


def draw_paired_estimates(axes, V, N, ST):
    """a to d: the four estimates, one line per individual from Illumina to Aviti, the median of each platform as a bar with its value,
    the Holm-corrected paired Wilcoxon result (genomescope_stats.tsv) as stars over the pair; the estimate as the title of each axis"""
    for ax, (e, lab), tag in zip(axes, EST, 'abcd'):
        il = np.array([float(V[('illumina', i)][e]) for i in range(1, N + 1)])
        av = np.array([float(V[('aviti', i)][e]) for i in range(1, N + 1)])
        for a, b in zip(il, av):
            ax.plot([0, 1], [a, b], color='#c5ccd3', lw=0.7, zorder=1)
        for x, v, (_, col) in ((0, il, PLAT['illumina']), (1, av, PLAT['aviti'])):
            ax.scatter(np.full(len(v), x), v, s=9, color=col, zorder=3, lw=0)
            ax.plot([x - 0.22, x + 0.22], [np.median(v)] * 2, color=col, lw=2.2, zorder=4)
            ax.text(x + (0.26 if x else -0.26), np.median(v), f'{np.median(v):.2f}', ha='left' if x else 'right', va='center',
                    fontsize=6.3, fontweight='bold', color=col)
        lo, hi = min(il.min(), av.min()), max(il.max(), av.max())
        pad = (hi - lo) * 0.12
        ax.set_ylim(lo - pad, hi + pad * 2.2)
        p = float(ST[e]['p_holm'])
        ax.plot([0, 0, 1, 1], [hi + pad * 1.2, hi + pad * 1.45, hi + pad * 1.45, hi + pad * 1.2], color=INK, lw=0.7)
        ax.text(0.5, hi + pad * 1.5, stars(p), ha='center', va='bottom', fontsize=7, color=INK)
        ax.set_xlim(-1.0, 1.6)          # room for the Illumina median value left of its bar, inside the axis line
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Illumina', 'Aviti'])
        ax.tick_params(axis='x', length=0)
        ax.set_title(lab, fontsize=7.5)
        ax.grid(axis='y', color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.text(-0.32, 1.06, tag, transform=ax.transAxes, fontsize=11, fontweight='bold', va='bottom')


def draw_het_vs_length(ax, V, N, COR):
    """the heterozygosity against the genome length estimate, one point per individual and platform; the Spearman rho of each platform
    (genomescope_correlations.tsv) in the key"""
    hand = []
    for plat in ('illumina', 'aviti'):
        name, col = PLAT[plat]
        h = [float(V[(plat, i)]['heterozygosity']) for i in range(1, N + 1)]
        g = [float(V[(plat, i)]['genome_length_gb']) for i in range(1, N + 1)]
        ax.scatter(g, h, s=14, color=col, lw=0)
        rho = next(float(r['spearman_rho']) for r in COR if r['platform'] == name and r['y'] == 'genome_length_gb')
        hand.append(Line2D([], [], color=col, marker='o', ls='', ms=4, label=f'{name} (Spearman rho {rho:.2f})'))
    ax.set_xlabel('Genome length estimate (Gb)')
    ax.set_ylabel('Heterozygosity (%)')
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(handles=hand, loc='upper right', fontsize=6.5)
