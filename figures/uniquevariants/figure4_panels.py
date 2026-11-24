#!/usr/bin/env python3
"""figure4_panels.py - the panels of Figure 4 and the drawing functions that figure4_combine.py and the supplementary scripts import:
a unique variants per individual by region (violins), b allele frequency in gnomAD genomes (numbers of alleles), c novel alleles by
region, d PPR-only variants (vs GATK) also called by DV. Panel e (Ti/Tv and het/hom by region) is still written as a single panel but is
not part of Figure 4; its function draws Supplementary Figure S18. Colours: Paul Tol's muted scheme, dark Illumina, light Aviti.
Inputs: the uv_*.tsv tables of uniquevariants_table.py (same folder). Outputs: figures/Figure4a-e (.pdf, .png at 600 dpi)
"""
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REGIONS = ['BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']   # the order of Figure 1
ROWS = ['genome'] + REGIONS
LAB = dict({r: r for r in REGIONS}, genome='Genome')
# regions by size of the counts, in the order of Figure 1 within a group
GROUPS = [['genome', 'BB', 'TRHP', 'LM', 'OD', 'AD'], ['15GC85', '25GC65', 'OMIM'], ['REFSEQ', 'ACMG']]
TYPES = [('ALL', 'All variants'), ('SNV', 'SNVs'), ('INDEL', 'Indels')]
PLATFORMS = [('illumina', 'Illumina'), ('aviti', 'Aviti')]
# the colours of Figure 1 (Paul Tol's muted scheme): (Illumina, Aviti); grey for the variants shared by two pipelines
SHADE = {'GATK': ('#882255', '#CC6677'), 'DV': ('#332288', '#88CCEE'), 'PPR': ('#117733', '#44AA99'), 'shared': ('#6E6E6E', '#BDBDBD')}
PANEL_TITLE = {'a': 'Pipeline-unique variants per individual', 'b': 'Pipeline-unique variants found in gnomAD genomes',
               'c': 'Novel variants by region',
               'd': 'PPR-only variants also called by DV', 'e': 'Ti/Tv and het/hom by region'}
D_YLABEL = 'PPR-only\nvariants vs\nGATK (%)'
PAIRS = [('GATK', 'PPR'), ('DV', 'PPR')]
HATCH = {'GATK': '', 'DV': '////'}            # bars of the pair DV vs PPR are hatched
INK, GRID = '#1f2933', '#e3e7eb'
EDGES = np.round(np.arange(-6.0, 0.0001, 0.25), 2)
CODING = ['missense', 'synonymous', 'splice', 'frameshift', 'in-frame', 'stop gained, stop lost, start lost']
CODING_LABEL = {'stop gained, stop lost, start lost': 'stop gained,\nstop or start lost'}
RNG = np.random.default_rng(17)


def style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'axes.titlesize': 7.5, 'axes.labelsize': 7,
                         'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'axes.edgecolor': '#9aa5b1', 'axes.labelcolor': INK,
                         'xtick.color': INK, 'ytick.color': INK, 'axes.spines.top': False, 'axes.spines.right': False,
                         'legend.frameon': False, 'axes.unicode_minus': False, 'hatch.linewidth': 0.5})


def read(name, folder=HERE):
    with open(os.path.join(folder, name), newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def load(folder=HERE):
    """C: (set, platform, region, type) -> 23 counts (individual 1-23); H: (set, platform, region, type, source) -> {bin_lo: alleles};
    Q: (set, platform, region, kind, group) -> alleles"""
    raw = defaultdict(dict)
    for r in read('uv_counts.tsv', folder):
        raw[(r['set'], r['platform'], r['region'], r['type'])][int(r['individual'])] = int(r['n'])
    C = {}
    for k, v in raw.items():
        if sorted(v) != list(range(1, 24)):
            raise SystemExit(f'{k}: {len(v)} individuals')
        C[k] = np.array([v[i] for i in range(1, 24)], dtype=float)
    H = defaultdict(dict)
    for r in read('uv_af_hist.tsv', folder):
        H[(r['set'], r['platform'], r['region'], r['type'], r['source'])][r['bin_lo']] = int(r['alleles'])
    Q = {(r['set'], r['platform'], r['region'], r['kind'], r['group']): int(r['alleles']) for r in read('uv_consequences.tsv', folder)}
    return C, H, Q


def set_label(a, b):
    return 'PPR only (vs GATK and DV)' if b == 'BOTH' else f'{a} only (vs {b})'


def kde(t, grid):
    """Gaussian kernel density on grid (Scott bandwidth), scaled to a maximum of 1"""
    sd = float(np.std(t, ddof=1)) if len(t) > 1 else 0.0
    bw = max(1.06 * sd * len(t) ** (-0.2), 0.01 * float(np.median(t)), 1.0)
    d = np.exp(-0.5 * ((grid[:, None] - t[None, :]) / bw) ** 2).sum(axis=1)
    return d / d.max(), bw


def half_violin(ax, x0, values, side, color, width=0.19):
    """one half of a violin (left = Illumina, right = Aviti) on a linear axis; median as a black line; individuals as small points"""
    t = np.asarray(values, dtype=float)
    _, bw = kde(t, t)
    grid = np.clip(np.linspace(t.min() - 1.2 * bw, t.max() + 1.2 * bw, 160), 0, None)
    d, _ = kde(t, grid)
    sgn = -1 if side == 'left' else 1
    ax.fill_betweenx(grid, x0, x0 + sgn * width * d, color=color, lw=0.4, ec='white', zorder=2)
    med = float(np.median(t))
    ax.plot([x0, x0 + sgn * width * np.interp(med, grid, d)], [med, med], color='black', lw=1.0, zorder=4)
    jit = RNG.uniform(0.12, 0.62, len(t)) * width * np.interp(t, grid, d)
    ax.scatter(x0 + sgn * jit, t, s=1.2, color=INK, alpha=0.45, lw=0, zorder=3)


def thousands(v, _):
    return f'{v / 1000:g}k' if v >= 1000 else f'{v:g}'


def draw_violins(axes, C, pair, typ='ALL', labels=True):
    """panel a for one pair: three axes (the groups of regions), two split violins per region (one per pipeline)"""
    a, b = pair
    for ax, rows in zip(axes, GROUPS):
        gx = np.arange(len(rows), dtype=float)
        for x, r in zip(gx, rows):
            for pipe, other, off in ((a, b, -0.21), (b, a, 0.21)):
                for (plat, _), side in zip(PLATFORMS, ('left', 'right')):
                    half_violin(ax, x + off, C[(f'{pipe}vs{other}', plat, r, typ)], side, SHADE[pipe][0 if plat == 'illumina' else 1])
        ax.set_xticks(gx)
        ax.set_xticklabels([LAB[r] for r in rows] if labels else [])
        ax.tick_params(axis='x', length=2)
        ax.set_xlim(-0.6, len(rows) - 0.4)
        ax.set_ylim(0, None)
        ax.grid(axis='y', color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.yaxis.set_major_formatter(FuncFormatter(thousands))
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4))


def legend_violins():
    return [Patch(color=SHADE[p][0], label=f'{p} only') for p in ('GATK', 'DV', 'PPR')] + [
        Patch(facecolor='#555555', label='left half: Illumina'), Patch(facecolor='#aaaaaa', label='right half: Aviti'),
        Line2D([], [], color='black', lw=1.0, label='median of the 23 individuals')]


AF_ABSENT = (-6.55, -6.3)     # with_absent: x of the "absent" bars of the two pipelines, left of the lowest AF bin (-5.5)


def draw_af(ax, H, pair, plat, source='gnomADg', region='genome', typ='ALL', ylabel=True, with_absent=False):
    """panel b: histograms of log10 AF of the unique alleles found in the source, as numbers of alleles; the share of the unique
    alleles found in the source is written in the panel. with_absent (not used in the figures): also one bar per pipeline for the
    alleles absent from the source, left of the histogram, labelled 'absent'"""
    a, b = pair
    notes = []
    for k, (pipe, other) in enumerate(((a, b), (b, a))):
        h = H[(f'{pipe}vs{other}', plat, region, typ, source)]
        present = np.array([h[f'{lo:.2f}'] for lo in EDGES[:-1]], dtype=float)
        c = SHADE[pipe][0]
        ax.stairs(present, EDGES, fill=True, color=c, alpha=0.35)
        ax.stairs(present, EDGES, color=c, lw=1.0)
        if with_absent:
            ax.bar(AF_ABSENT[k], h['absent'], width=0.25, color=c, alpha=0.35, edgecolor=c, lw=1.0)
        notes.append((c, f'{pipe} only {100 * present.sum() / (present.sum() + h["absent"]):.0f}%'))
    ax.yaxis.set_major_formatter(FuncFormatter(thousands))
    if with_absent:
        ax.set_xlim(-6.8, 0)
        ax.axvline(-5.85, color='#9aa5b1', lw=0.6, ls=(0, (2, 2)))
        ax.set_xticks([sum(AF_ABSENT) / 2, -4, -2, 0])
        ax.set_xticklabels(['absent', '0.0001', '0.01', '1'])
    else:
        ax.set_xlim(-5.5, 0)
        ax.set_xticks([-4, -2, 0])
        ax.set_xticklabels(['0.0001', '0.01', '1'])
    ax.set_xticks(range(-5, 1), minor=True)
    ax.grid(axis='y', color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    # the notes go to the right when the absent bars stand on the left
    nx, ha = (0.97, 'right') if with_absent else (0.03, 'left')
    ax.text(nx, 0.97, 'found in gnomAD:' if source == 'gnomADg' else 'found in 1000 Genomes:', transform=ax.transAxes, va='top',
            ha=ha, fontsize=plt.rcParams['xtick.labelsize'], color=INK)
    for k, (c, t) in enumerate(notes):
        ax.text(nx, 0.86 - 0.11 * k, t, transform=ax.transAxes, va='top', ha=ha, fontsize=plt.rcParams['xtick.labelsize'], color=c,
                fontweight='bold')
    if ylabel:
        found = 'found in gnomAD' if source == 'gnomADg' else 'found in 1000 Genomes'
        ax.set_ylabel('Pipeline-unique variants\n' + ('in the cohort' if with_absent else found))


def draw_consequences(ax, Q, region='genome', values=True):
    """panel c: coding and splice consequences, distinct alleles of the cohort, four sides x two platforms per class"""
    sides = [(a, b) for pair in PAIRS for a, b in (pair, pair[::-1])]
    n_bars = len(sides) * 2
    wbar = 0.84 / n_bars
    top = 0
    for g, cl in enumerate(CODING):
        k = 0
        for a, b in sides:
            hatch = HATCH['DV'] if 'DV' in (a, b) else ''
            for pi, (plat, _) in enumerate(PLATFORMS):
                v = Q[(f'{a}vs{b}', plat, region, 'coding', cl)]
                top = max(top, v)
                x = g - 0.42 + wbar * (k + 0.5)
                ax.bar(x, v, width=wbar * 0.9, color=SHADE[a][pi], hatch=hatch, lw=0.4,
                       edgecolor=(SHADE[a][0] if pi else 'white') if hatch else 'white')
                if values:   # the value above each bar
                    ax.text(x, v + top * 0.012, f'{v:,}', rotation=90, ha='center', va='bottom', fontsize=6.0, fontweight='bold',
                            color=INK)
                k += 1
    ax.set_xticks(range(len(CODING)))
    ax.set_xticklabels([CODING_LABEL.get(c, c) for c in CODING])
    ax.tick_params(axis='x', length=0)
    ax.set_ylabel('Pipeline-unique variants\nin the cohort')
    ax.set_ylim(0, top * 1.4)   # room for the vertical values above the highest bar
    ax.set_xlim(-0.5, len(CODING) - 0.5)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:,.0f}'))
    ax.grid(axis='y', color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def legend_consequences():
    h = []
    for a, b in [(a, b) for pair in PAIRS for a, b in (pair, pair[::-1])]:
        hatch = HATCH['DV'] if 'DV' in (a, b) else ''
        h.append(Patch(facecolor=SHADE[a][0], edgecolor='white', hatch=hatch, label=f'{set_label(a, b)}, Illumina'))
        h.append(Patch(facecolor=SHADE[a][1], edgecolor=SHADE[a][0] if hatch else 'white', hatch=hatch, label=f'{set_label(a, b)}, Aviti'))
    return h


def split_shares(C, region, typ, plat):
    """per individual: the share of the PPR-only variants (vs GATK) that DV also misses (PPR only vs GATK and DV)"""
    g, both = C[('PPRvsGATK', plat, region, typ)], C[('PPRvsBOTH', plat, region, typ)]
    if np.any(both > g):
        raise SystemExit(f'{region} {typ} {plat}: PPR only vs GATK and DV above PPR only vs GATK')
    return both[g > 0] / g[g > 0]


def region_axis(ax, rows, labels=True):
    """the shared x axis of panels c, d and e: one group per region, in the order of Figure 1"""
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([('Autosomes' if r == 'autosomes' else LAB[r]) for r in rows])
    ax.tick_params(axis='x', length=0, labelbottom=labels)
    ax.set_xlim(-0.5, len(rows) - 0.5)
    ax.grid(axis='y', color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def draw_caller_graph(ax, C, typ, labels=True):
    """panel d: by region, the PPR-only variants (vs GATK) split into those BWA-MEM + DV also calls (lower part) and those only PPR
    calls (upper part); median of the 23 individuals; black line: interquartile range of the split; left bar Illumina, right bar Aviti"""
    for g, r in enumerate(ROWS):
        for pi, ((plat, _), off) in enumerate(zip(PLATFORMS, (-0.19, 0.19))):
            sh = split_shares(C, r, typ, plat)
            graph = 100 * float(np.median(sh))
            x = g + off
            ax.bar(x, 100 - graph, width=0.34, color=SHADE['DV'][pi], edgecolor='white', lw=0.3)
            ax.bar(x, graph, bottom=100 - graph, width=0.34, color=SHADE['PPR'][pi], edgecolor='white', lw=0.3)
            ax.vlines(x, 100 - 100 * np.percentile(sh, 75), 100 - 100 * np.percentile(sh, 25), color=INK, lw=1.0)
    region_axis(ax, ROWS, labels)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 50, 100])
    ax.set_ylabel(D_YLABEL)


def legend_caller_graph():
    """one row: the two parts in the Illumina shade, the two shades, the range"""
    return [Patch(color=SHADE['DV'][0], label='also called by BWA-MEM + DV'), Patch(color=SHADE['PPR'][0], label='called by PPR only'),
            Patch(facecolor='#555555', label='dark: Illumina'), Patch(facecolor='#bbbbbb', label='light: Aviti'),
            Line2D([], [], color=INK, lw=1.0, label='interquartile range')]


def load_r27(folder=HERE):
    """T: (pair, set or 'shared', platform, region) -> (Ti/Tv, het/hom) for the pairs GATK vs PPR and DV vs PPR (shared: the mean of
    the two sides, which differ only in the way a few alleles are written; nan where a ratio is NA); N: (set, platform, region) ->
    number of novel alleles of the cohort"""
    T = {}
    for r in read('uv_titv_hethom.tsv', folder):
        if r['pair'] in ('GATK_PPR', 'DV_PPR') and r['class'] in ('shared', 'unique'):
            key = (r['pair'], 'shared' if r['class'] == 'shared' else r['set'], r['platform'], r['region'])
            T.setdefault(key, []).append(tuple(float('nan') if r[c] in ('', 'NA') else float(r[c]) for c in ('titv', 'hethom')))
    T = {k: tuple(float(np.mean([v[i] for v in vs])) for i in (0, 1)) for k, vs in T.items()}
    N = {(r['set'], r['platform'], r['region']): int(r['novel']) for r in read('uv_novel.tsv', folder)}
    return T, N


NOVEL_SETS = [('GATK', 'PPR'), ('PPR', 'GATK'), ('DV', 'PPR'), ('PPR', 'DV')]
# the novel alleles were counted on the autosomes only, in every region too; regions in the groups of panel a (by size)
NOVEL_GROUPS = [['autosomes', 'BB', 'TRHP', 'LM', 'OD', 'AD'], ['15GC85', '25GC65', 'OMIM'], ['REFSEQ', 'ACMG']]
MARK = {'GATK_PPR': 'o', 'DV_PPR': 's'}          # circles: the pair GATK vs PPR; squares: the pair DV vs PPR


def draw_novel(axes, N):
    """panel c: the number of novel alleles (absent from dbSNP, gnomAD genomes and 1000 Genomes) among the unique alleles of the cohort,
    autosomal alleles, by region, as lollipops, the regions in three groups by size, each with its own axis; circles GATK vs PPR,
    squares DV vs PPR; dark Illumina, light Aviti"""
    w = 0.80 / (len(NOVEL_SETS) * len(PLATFORMS) + 0.5)
    for ax, rows in zip(axes, NOVEL_GROUPS):
        for g, reg in enumerate(rows):
            k = 0
            for a, b in NOVEL_SETS:
                mark = MARK['DV_PPR'] if 'DV' in (a, b) else MARK['GATK_PPR']
                for pi, (plat, _) in enumerate(PLATFORMS):
                    x = g - 0.40 + w * (k + 0.5 + (0.5 if k >= 4 else 0.0))
                    v = N[(f'{a}vs{b}', plat, reg)]
                    ax.plot([x, x], [0, v], color=SHADE[a][pi], lw=1.0, zorder=2)
                    ax.scatter([x], [v], s=20, marker=mark, color=SHADE[a][pi], edgecolor='white', lw=0.3, zorder=3)
                    k += 1
        ax.set_xticks(range(len(rows)))
        ax.set_xticklabels([('Autosomes' if r == 'autosomes' else LAB[r]) for r in rows])
        ax.tick_params(axis='x', length=0)
        ax.set_xlim(-0.6, len(rows) - 0.4)
        ax.set_ylim(0, None)
        ax.yaxis.set_major_formatter(FuncFormatter(thousands))
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4))
        ax.grid(axis='y', color=GRID, lw=0.6)
        ax.set_axisbelow(True)
    axes[0].set_ylabel('Novel variants\nin the cohort')


def legend_novel():
    """one row: the four sets in the Illumina shade (circles GATK vs PPR, squares DV vs PPR) and the two shades"""
    h = [Line2D([], [], marker=MARK['DV_PPR'] if 'DV' in (a, b) else MARK['GATK_PPR'], ls='', ms=4.5, color=SHADE[a][0],
                label=set_label(a, b)) for a, b in NOVEL_SETS]
    return h + [Line2D([], [], marker='o', ls='', ms=4.5, color='#555555', label='dark: Illumina'),
                Line2D([], [], marker='o', ls='', ms=4.5, color='#bbbbbb', label='light: Aviti')]


E_PAIRS = [('GATK_PPR', 'GATK', 'PPR'), ('DV_PPR', 'DV', 'PPR')]


def draw_titv_hethom(axes, T, labels=True):
    """panel e: Ti/Tv of the SNVs (first axis) and het/hom (second axis, log) by region, sums over the 23 individuals, for the pairs
    GATK vs PPR (left half of a region, circles) and DV vs PPR (right half, squares): the variants found by one pipeline of the pair as
    dots, those shared by both as a grey line across the half; dark Illumina, light Aviti"""
    offs = (-0.12, -0.04, 0.04, 0.12)
    for ax, what in zip(axes, ('titv', 'hethom')):
        i = 0 if what == 'titv' else 1
        for g, reg in enumerate(ROWS):
            for j, (pair, a, b) in enumerate(E_PAIRS):
                x0 = g + (-0.21 if j == 0 else 0.21)
                for pi, (plat, _) in enumerate(PLATFORMS):
                    v = T[(pair, 'shared', plat, reg)][i]
                    ax.plot([x0 - 0.17, x0 + 0.17], [v, v], color=SHADE['shared'][pi], lw=1.4, zorder=1, solid_capstyle='butt')
                for k, s in enumerate((a, b)):
                    for pi, (plat, _) in enumerate(PLATFORMS):
                        ax.scatter([x0 + offs[2 * k + pi]], [T[(pair, s, plat, reg)][i]], s=13, marker=MARK[pair],
                                   color=SHADE[s][pi], edgecolor='white', lw=0.3, zorder=3)
        region_axis(ax, ROWS, labels and what == 'hethom')
        ax.grid(axis='x', color=GRID, lw=0.4)
        if what == 'hethom':
            ax.set_yscale('log')
            ax.set_ylim(0.25, 32)
            ax.set_yticks([0.3, 1, 3, 10, 30])
            ax.set_yticklabels(['0.3', '1', '3', '10', '30'])
            ax.minorticks_off()
            ax.axhline(1, color='#9aa5b1', lw=0.6, ls=(0, (3, 2)))
        else:
            ax.set_ylim(0.6, 3.2)
            ax.set_yticks([1.0, 1.5, 2.0, 2.5, 3.0])
        ax.set_ylabel({'titv': 'Ti/Tv (SNVs)', 'hethom': 'het/hom'}[what])


def legend_r27():
    """one row: the shared line, the four sets in the Illumina shade (circles GATK vs PPR, squares DV vs PPR) and the two shades"""
    return [Line2D([], [], color=SHADE['shared'][0], lw=1.4, label='shared by the pair')] + legend_novel()


def panel_tags(fig, items, dy=0.012):
    """the tags and short titles of the panels of a supplementary figure, items = [(axes, tag, title)]: above the top of each axes, at
    the left edge of everything already drawn (call it after drawing, before supp_title)"""
    fig.canvas.draw()
    x0 = fig.get_tightbbox(fig.canvas.get_renderer()).x0 / fig.get_figwidth()
    for ax, tag, title in items:
        y = ax.get_position().y1 + dy
        fig.text(x0, y, tag, fontsize=11, fontweight='bold', va='bottom')
        fig.text(x0 + 0.035, y + 0.001, title, fontsize=8, fontweight='bold', va='bottom')


def supp_title(fig, text):
    """the title of a supplementary figure: top left, above everything already drawn (call it just before save)"""
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())   # inches
    fig.text(bb.x0 / fig.get_figwidth(), (bb.y1 + 0.08) / fig.get_figheight(), text, ha='left', va='bottom', fontsize=8.5,
             fontweight='bold', color=INK)


def save(fig, stem):
    os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True)
    for ext, kw in (('pdf', {}), ('png', {'dpi': 600})):
        fig.savefig(os.path.join(HERE, 'figures', f'{stem}.{ext}'), bbox_inches='tight', facecolor='white', **kw)
    plt.close(fig)
    print('written', stem)


if __name__ == '__main__':
    style()
    C, H, Q = load()
    fig = plt.figure(figsize=(7.1, 3.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[6, 3, 2], wspace=0.25, hspace=0.35)
    for i, pair in enumerate(PAIRS):
        axes = [fig.add_subplot(gs[i, j]) for j in range(3)]
        draw_violins(axes, C, pair)
        axes[0].text(0.0, 1.03, f'{pair[0]} vs {pair[1]}', transform=axes[0].transAxes, fontweight='bold', fontsize=7.5, va='bottom')
    fig.supylabel('Pipeline-unique variants per individual', fontsize=7, x=0.05)
    fig.legend(handles=legend_violins(), loc='upper center', bbox_to_anchor=(0.5, 0.03), ncol=6, fontsize=6)
    save(fig, 'Figure4a_violins')
    fig, axes = plt.subplots(1, 4, figsize=(7.1, 1.6), sharex=True)
    for j, (pair, (plat, pl)) in enumerate([(p, q) for p in PAIRS for q in PLATFORMS]):
        draw_af(axes[j], H, pair, plat, ylabel=(j == 0))
        axes[j].set_title(f'{pair[0]} vs {pair[1]}, {pl}', fontweight='bold')
    fig.supxlabel('Allele frequency in gnomAD genomes', fontsize=7, y=-0.12)
    fig.subplots_adjust(wspace=0.3)
    save(fig, 'Figure4b_gnomad')
    T, N = load_r27()
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 1.7), gridspec_kw={'width_ratios': [6, 3, 2], 'wspace': 0.22})
    draw_novel(axes, N)
    axes[0].set_title(PANEL_TITLE['c'], fontweight='bold', fontsize=7.5, loc='left')
    fig.legend(handles=legend_novel(), loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=6, fontsize=5.8)
    save(fig, 'Figure4c_novel')
    fig, ax = plt.subplots(figsize=(7.1, 1.3))
    draw_caller_graph(ax, C, 'ALL')
    ax.set_title(PANEL_TITLE['d'], fontweight='bold', fontsize=7.5, loc='left')
    fig.legend(handles=legend_caller_graph(), loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=5, fontsize=5.8)
    save(fig, 'Figure4d_caller_graph')
    fig, axes = plt.subplots(2, 1, figsize=(7.1, 2.8), sharex=True)
    draw_titv_hethom(axes, T)
    axes[0].set_title(PANEL_TITLE['e'], fontweight='bold', fontsize=7.5, loc='left')
    fig.subplots_adjust(hspace=0.12)
    fig.legend(handles=legend_r27(), loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=7, fontsize=5.8)
    save(fig, 'Figure4e_titv_hethom')
    print('FIGURE4_PANELS_DONE')
