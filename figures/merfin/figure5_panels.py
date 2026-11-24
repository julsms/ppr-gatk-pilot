#!/usr/bin/env python3
"""figure5_panels.py - drawing functions of Figure 5 and of the Merfin supplementary figures (imported), with the tables of
merfin_table.py from this folder. Colours: Paul Tol's muted scheme as in Figure 1 (dark shade Illumina, light shade Aviti); GIAB
classes in Tol grey, purple and sand; olive for Figure 5b and 5d; three context colours for Figure 5c.
"""
import csv
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REGIONS = ['BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']   # the order of Figure 1
ROWS = ['genome'] + REGIONS
LAB = dict({r: r for r in REGIONS}, genome='Genome')
PLATS = [('illumina', 'Illumina'), ('aviti', 'Aviti')]
# the colours of Figure 1 (Paul Tol's muted scheme): (Illumina, Aviti)
SHADE = {'GATK': ('#882255', '#CC6677'), 'DV': ('#332288', '#88CCEE'), 'PPR': ('#117733', '#44AA99')}
FULL = [('GATK_ALL', 'GATK', 'GATK'), ('DV_ALL', 'DV', 'DV'), ('PPR_ALL', 'PPR', 'PPR')]
SPEC = [('GATKvsPPR', 'GATK only\n(vs PPR)', 'GATK'), ('PPRvsGATK', 'PPR only\n(vs GATK)', 'PPR'),
        ('DVvsPPR', 'DV only\n(vs PPR)', 'DV'), ('PPRvsDV', 'PPR only\n(vs DV)', 'PPR'),
        ('PPRvsBOTH', 'PPR only\n(vs GATK, DV)', 'PPR')]
SETS = FULL + SPEC
LABEL = {s: lab for s, lab, _ in SETS}
PIPE = {s: p for s, _, p in SETS}
HATCH = {'DVvsPPR': '////', 'PPRvsDV': '////', 'PPRvsBOTH': 'xxxx'}
INK, GRID, MUTED = '#1f2933', '#e3e7eb', '#7b8794'
CAT = {'TP': '#BBBBBB', 'FP': '#AA4499', 'FP_gt': '#DDCC77'}          # GIAB classes (Tol grey, purple, sand)
# k-mer context: the classes that k-mers cannot place (darker = more copies), then single copy
CTX = [('kstar', 'Excess k-mer copies', '#736600'), ('repeat', 'Repeated in GRCh38, no excess copies', '#C8B43B'),
       ('single', 'Single-copy context', '#CED6EF')]
OLIVE, OLIVE_LIGHT = '#736600', '#C4BE92'      # Figure 5: b and d; a: Illumina, Aviti
OPEN = '#f5f7fa'                               # open bars: GIAB precision
SHARE_GREY = '#616e7c'                         # Figure 5d, percentage inside the benchmark: contrast 5.2:1 on white
P6SETS = [('PPRvsGATK', 'PPR only (vs GATK)'), ('GATKvsPPR', 'GATK only (vs PPR)')]
PANEL_TITLE = {'a': 'K-mer support in the cohort', 'b': 'K-mer support and GIAB precision, HG002',
               'c': 'K-mer context of the supported HG002 records', 'd': 'K-mer support and GIAB precision by region, HG002'}


def style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.titlesize': 8.5, 'axes.labelsize': 7.5,
                         'xtick.labelsize': 6.8, 'ytick.labelsize': 7, 'axes.edgecolor': '#9aa5b1', 'axes.labelcolor': INK,
                         'xtick.color': INK, 'ytick.color': INK, 'axes.spines.top': False, 'axes.spines.right': False,
                         'legend.frameon': False, 'axes.unicode_minus': False, 'hatch.linewidth': 0.5})


def read(name, folder=HERE):
    with open(os.path.join(folder, name), newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def num(v):
    return float('nan') if v in ('', 'NA') else float(v)


def load(folder=HERE):
    """R: (platform, set, region, type) -> Merfin rate (%) of the 23 individuals; NV: (platform, set, region) -> (n, validated) of the
    23 individuals, all variants; HG: (set, platform) -> the HG002 numbers; HR: (platform, set, region) -> the HG002 numbers by region"""
    raw = {}
    for r in read('mf_rates.tsv', folder):
        raw.setdefault((r['platform'], r['set'], r['region'], r['type']), {})[int(r['individual'])] = (int(r['n']), int(r['validated']))
    R, NV = {}, {}
    for k, v in raw.items():
        if sorted(v) != list(range(1, 24)):
            raise SystemExit(f'{k}: {len(v)} individuals')
        n = np.array([v[i][0] for i in range(1, 24)], dtype=float)
        val = np.array([v[i][1] for i in range(1, 24)], dtype=float)
        with np.errstate(invalid='ignore', divide='ignore'):
            R[k] = np.where(n > 0, 100 * val / n, np.nan)
        if k[3] == 'ALL':
            NV[k[:3]] = (n, val)
    HG = {(r['set'], r['platform']): {c: num(r[c]) for c in r if c not in ('set', 'platform')} for r in read('mf_hg002.tsv', folder)}
    HR = {(r['platform'], r['set'], r['region']): {c: num(r[c]) for c in ('records', 'evaluable', 'precision', 'merfin', 'inside')}
          for r in read('mf_hg002_regions.tsv', folder)}
    return R, NV, HG, HR


def save(fig, stem, folder=HERE):
    os.makedirs(os.path.join(folder, 'figures'), exist_ok=True)
    for ext, kw in (('pdf', {}), ('png', {'dpi': 600})):
        fig.savefig(os.path.join(folder, 'figures', f'{stem}.{ext}'), bbox_inches='tight', facecolor='white', **kw)
    plt.close(fig)
    print('written', stem)


def supp_title(fig, text):
    """the title of a supplementary figure: top left, above everything already drawn (call it just before save)"""
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    fig.text(bb.x0 / fig.get_figwidth(), (bb.y1 + 0.08) / fig.get_figheight(), text, ha='left', va='bottom', fontsize=8.5,
             fontweight='bold', color=INK)


def bar(ax, x, h, sid, pi, width=0.72, **kw):
    """a bar in the colour of the pipeline that made the set; the DV vs PPR pair hatched, PPR only vs both cross-hatched"""
    hatch = HATCH.get(sid, '')
    edge = (SHADE[PIPE[sid]][0] if pi else 'white') if hatch else 'white'
    ax.bar(x, h, width=width, color=SHADE[PIPE[sid]][pi], hatch=hatch, edgecolor=edge, lw=0.5, **kw)


def xpos():
    """x of the 8 sets: three full sets, a gap, five sets of one pipeline"""
    return {s: (i if i < 3 else i + 0.6) for i, (s, _, _) in enumerate(SETS)}


def set_axis(ax, fs, labels=True):
    X = xpos()
    ax.set_xticks([X[s] for s in LABEL])
    ax.set_xticklabels([LABEL[s] for s in LABEL] if labels else [], rotation=90, fontsize=fs, linespacing=1.0)
    ax.tick_params(axis='x', length=0)
    ax.set_xlim(-0.5, max(X.values()) + 0.5)      # room for the two-line names at 7 pt, side by side


def legend_sets(pi=0):
    """the colours and hatches of the sets of one pipeline, in one shade"""
    return [Patch(facecolor=SHADE[PIPE[s]][pi], edgecolor=(SHADE[PIPE[s]][0] if pi else 'white') if HATCH.get(s) else 'white',
                  hatch=HATCH.get(s, ''), label=LABEL[s].replace('\n', ' ')) for s, _, _ in SPEC]


# ---------------------------------------------------------------- Figure 5
def draw_support_bars(ax, R, fs=7.0):
    """a: the mean percentage of the records with k-mer support in the cohort, genome, all variants; one row per set, two horizontal
    bars (Illumina dark olive, Aviti light olive; hatches as in b), the value at the end of each bar; the axis starts at zero;
    returns the legend handles"""
    sets = list(LABEL)
    h = 0.38
    for i, s in enumerate(sets):
        for pi, (plat, _) in enumerate(PLATS):
            v = float(np.nanmean(R[(plat, s, 'genome', 'ALL')]))
            y = i + (pi - 0.5) * h
            hatch = HATCH.get(s, '')
            ax.barh(y, v, height=h * 0.92, color=OLIVE_LIGHT if pi else OLIVE, hatch=hatch, lw=1.0 if hatch and pi else 0.5,
                    edgecolor=(OLIVE if pi else 'white') if hatch else 'white')
            ax.text(v + 1.5, y, f'{v:.1f}', ha='left', va='center', fontsize=fs, color=INK)   # outside: hatches stay legible
    ax.set_yticks(range(len(sets)))
    ax.set_yticklabels([LABEL[s].replace('\n', ' ') for s in sets], fontsize=fs)
    ax.invert_yaxis()
    ax.set_xlim(0, 118)
    ax.spines['bottom'].set_bounds(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.tick_params(axis='x', labelsize=fs)
    ax.tick_params(axis='y', length=0)
    ax.set_xlabel('Records with k-mer support (%)', fontsize=fs + 1)
    ax.grid(axis='x', color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    return [Patch(facecolor=OLIVE, edgecolor='white', label='Illumina'),
            Patch(facecolor=OLIVE_LIGHT, edgecolor='white', label='Aviti')]


MERFIN_GIAB = [Patch(facecolor=MUTED, label='With k-mer support (Merfin)'),
               Patch(facecolor=OPEN, edgecolor=INK, label='True positive (GIAB precision)')]


def draw_merfin_vs_giab(axes, HG, fs=6.3, vertical=False, col=None):
    """b: HG002, records inside the benchmark: percentage with k-mer support (filled bar) and GIAB precision (open bar) of the same
    records; the difference in points above each pair. col: one colour for every set (Figure 5; the sets are named on the axis and
    the DV pair keeps its hatching); None: the pipeline colours"""
    X = xpos()
    w = 0.38
    for ax, (plat, pl), pi in zip(axes, PLATS, (0, 1)):
        for s in LABEL:
            h = HG[(s, plat)]
            if col:
                ax.bar(X[s] - w / 2, h['merfin'], width=w, color=col, hatch=HATCH.get(s, ''), edgecolor='white', lw=0.5)
            else:
                bar(ax, X[s] - w / 2, h['merfin'], s, pi, width=w)
            ax.bar(X[s] + w / 2, h['precision'], width=w, color=OPEN, edgecolor=col or SHADE[PIPE[s]][0], lw=1.0)
            ax.text(X[s], max(h['merfin'], h['precision']) + 1.5, f'+{h["gap"]:.1f}', ha='center', va='bottom', fontsize=fs,
                    fontweight='bold', color=INK, rotation=90 if vertical else 0)
        set_axis(ax, fs)
        ax.set_ylim(0, 122 if vertical else 110)
        ax.set_yticks(range(0, 101, 20))
        ax.grid(axis='y', color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.set_title(f'HG002, {pl}', fontsize=fs + 1)
    axes[0].set_ylabel('Records inside the GIAB benchmark (%)')
    if col:
        return [Patch(facecolor=col, label='With k-mer support (Merfin)'),
                Patch(facecolor=OPEN, edgecolor=col, label='True positive (GIAB precision)')]
    return MERFIN_GIAB


def ypos():
    """rows of panel c, top to bottom: the three full sets, a gap, the five sets of one pipeline"""
    return {s: (i if i < 3 else i + 0.6) for i, (s, _, _) in enumerate(SETS)}


def not_single(h, g):
    """percentage of the supported records of a GIAB group (TP; FP = FP + FP_gt) outside a single-copy context"""
    return h[f'ctx_{g}_kstar'] + h[f'ctx_{g}_repeat']


def fits(ax, text, x0, x1, pad=3.0):
    """True when the drawn text is narrower than the data interval x0-x1 (with a few pixels to spare)"""
    width = text.get_window_extent(renderer=ax.figure.canvas.get_renderer()).width
    return width + pad <= abs(ax.transData.transform((x1, 0))[0] - ax.transData.transform((x0, 0))[0])


def draw_context(axes, HG, fs=6.3):
    """c: HG002, the supported records that GIAB confirms (TP, left of the centre line) or rejects (FP, right), by the k-mer context
    of the allele, the classes from the centre out; on each bar the percentage outside a single-copy context and, on the light part,
    the single-copy percentage, in bold; a single-copy value that does not fit inside its part is written just after the end of the
    bar, and a value outside a single-copy context that would run into the dark part starts just outside the dark part"""
    Y = ypos()
    for ax, (plat, pl) in zip(axes, PLATS):
        ax.set_xlim(-124, 124)                    # room for a 7-pt value written after the end of a bar
        for s in LABEL:
            h = HG[(s, plat)]
            for side, g in ((-1, 'TP'), (1, 'FP')):
                base = 0.0
                for c, _, col in CTX:
                    v = h[f'ctx_{g}_{c}']
                    ax.barh(Y[s], side * v, left=side * base, height=0.68, color=col, edgecolor='white', lw=0.4)
                    base += v
                ns, single = not_single(h, g), h[f'ctx_{g}_single']
                t = ax.text(side * (ns - 1.5), Y[s], f'{ns:.1f}', ha='right' if side > 0 else 'left', va='center', fontsize=fs,
                            fontweight='bold', color=INK)
                if not fits(ax, t, side * h[f'ctx_{g}_kstar'], side * (ns - 1.5)):
                    t.set_position((side * (h[f'ctx_{g}_kstar'] + 0.8), Y[s]))
                    t.set_ha('left' if side > 0 else 'right')
                t = ax.text(side * (ns + single / 2), Y[s], f'{single:.1f}', ha='center', va='center', fontsize=fs,
                            fontweight='bold', color=INK)
                if not fits(ax, t, side * ns, side * 100):
                    t.set_position((side * 101.5, Y[s]))
                    t.set_ha('left' if side > 0 else 'right')
        ax.axvline(0, color=INK, lw=1.0, zorder=3)
        ax.set_xticks([-100, -50, 0, 50, 100])
        ax.set_xticklabels(['100', '50', '0', '50', '100'], fontsize=fs)
        ax.set_yticks([Y[s] for s in LABEL])
        ax.set_yticklabels([LABEL[s].replace('\n', ' ') for s in LABEL], fontsize=fs)
        ax.set_ylim(max(Y.values()) + 0.6, -0.6)
        ax.tick_params(axis='y', length=0)
        ax.spines['left'].set_visible(False)
        ax.set_title(f'HG002, {pl}', fontsize=fs + 1, pad=11)
        tr = ax.get_xaxis_transform()
        ax.text(-50, 1.01, 'GIAB true positives', transform=tr, ha='center', va='bottom', fontsize=fs)
        ax.text(50, 1.01, 'GIAB false positives', transform=tr, ha='center', va='bottom', fontsize=fs)
        ax.set_xlabel('Records with k-mer support (%)', fontsize=fs + 1)
    return [Patch(facecolor=col, edgecolor='#9aa5b1', label=lab) for _, lab, col in CTX]


LINES = [Line2D([], [], color=INK, marker='o', ms=3.6, lw=1.4, label='With k-mer support (Merfin)'),
         Line2D([], [], color=INK, marker='o', ms=3.6, mfc='white', ls='--', lw=1.2, label='True positive (GIAB precision)')]


def draw_regions_lines(axes, HR, plat, side_by_side=False, fs=6.8, col=None):
    """d: HG002, by region: k-mer support (solid) and GIAB precision (dashed) of the same records, the gap shaded, for PPR only (vs
    GATK) and GATK only (vs PPR); the percentage of the records inside the benchmark in grey under the points. col: one colour for
    both sets (Figure 5); None: the colour of the pipeline that made the set (supplementary figure, Aviti). side_by_side (Figure 5d):
    the differences in one row above the data, each over its region, without a + sign (the script stops if one is negative), and the
    percentages inside the benchmark without the % sign; otherwise each difference just above the higher point"""
    pi = 0 if plat == 'illumina' else 1
    x = np.arange(len(ROWS))
    one = col
    for ax, (s, title) in zip(axes, P6SETS):
        col = one or SHADE[PIPE[s]][pi]
        r = [HR[(plat, s, st)] for st in ROWS]
        merf = np.array([v['merfin'] for v in r])
        prec = np.array([v['precision'] for v in r])
        ax.fill_between(x, prec, merf, color=col, alpha=0.14, lw=0)
        ax.plot(x, merf, '-o', color=col, lw=1.3 if side_by_side else 1.4, ms=3.6)
        ax.plot(x, prec, '--o', color=col, lw=1.1 if side_by_side else 1.2, ms=3.6, mfc='white')
        ax.set_ylim(0, 117 if side_by_side else 115)
        if side_by_side and np.any(merf < prec):
            raise SystemExit(f'{s}: a negative difference, which the labels without a sign cannot show')
        for i in range(len(ROWS)):
            y = 104 if side_by_side else max(merf[i], prec[i]) + 3
            ax.text(x[i], y, f'{merf[i] - prec[i]:.1f}' if side_by_side else f'+{merf[i] - prec[i]:.1f}', ha='center', va='bottom',
                    fontsize=fs if side_by_side else 5.8, fontweight='bold', color=INK)
            if side_by_side:
                ax.text(x[i], -3, f'{r[i]["inside"]:.1f}', ha='center', va='top', fontsize=fs, color=SHARE_GREY, clip_on=False)
            else:
                ax.text(x[i], -3, f'{r[i]["inside"]:.1f}%', ha='center', va='top', fontsize=5.6, color=MUTED, clip_on=False)
        ax.set_yticks(range(0, 101, 20))
        ax.grid(axis='y', color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.set_title(f'HG002, {dict(PLATS)[plat]}: {title}', fontsize=fs + 1 if side_by_side else 8, loc='left')
        ax.set_xticks(x)
        ax.set_xlim(-0.65 if side_by_side else -0.5, len(ROWS) - (0.35 if side_by_side else 0.5))   # d: the end labels inside
        ax.tick_params(axis='x', length=0, pad=12 if side_by_side else 11)
        ax.set_xticklabels([LAB[r] for r in ROWS], rotation=90 if side_by_side else 0, fontsize=fs if side_by_side else fs - 0.6)
    axes[0].set_ylabel('Records inside the\nGIAB benchmark (%)')
    return LINES


# ---------------------------------------------------------------- supplementary figures
def half_violin(ax, x, v, side, color, width=0.36):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    if len(v) >= 2 and v.std() > 0:
        bw = 1.06 * v.std(ddof=1) * len(v) ** (-0.2)
        g = np.linspace(v.min(), v.max(), 120)
        d = np.exp(-0.5 * ((g[:, None] - v[None, :]) / bw) ** 2).sum(1)
        d = width * d / d.max()
        if side < 0:
            ax.fill_betweenx(g, x - d, x, color=color, lw=0)
        else:
            ax.fill_betweenx(g, x, x + d, color=color, lw=0)
    med = np.median(v)
    ax.plot([x, x + side * width * 0.9], [med, med], color=INK, lw=0.9)
    jit = np.linspace(0.06, 0.2, len(v)) * side
    ax.scatter(x + jit, v, s=2.5, color=INK, alpha=0.55, lw=0, zorder=3)


def draw_violins(axes, R):
    """the Merfin rate of each individual: all calls (left axis) and the sets of one pipeline (right axis); left half Illumina,
    right half Aviti"""
    for ax, group, step in ((axes[0], FULL, 1.0), (axes[1], SPEC, 10.0)):
        low = 100.0
        for i, (s, lab, pipe) in enumerate(group):
            for (plat, _), side, pi in zip(PLATS, (-1, 1), (0, 1)):
                v = R[(plat, s, 'genome', 'ALL')]
                low = min(low, float(np.nanmin(v)))
                half_violin(ax, i, v, side, SHADE[pipe][pi])
        ax.set_xticks(range(len(group)))
        ax.set_xticklabels([lab for _, lab, _ in group], fontsize=6.5)
        ax.tick_params(axis='x', length=0)
        ax.set_ylim(step * np.floor(low / step) - (0.5 if step == 1.0 else 0.0), 100)   # every individual inside the axis
        ax.grid(axis='y', color=GRID, lw=0.6)
        ax.set_axisbelow(True)
    axes[0].set_ylabel('Records with k-mer support (%)')
    axes[0].set_title('All calls', fontsize=8)
    axes[1].set_title('Variants found by one pipeline', fontsize=8)
    return [Patch(color=SHADE[p][0], label=p) for p in ('GATK', 'DV', 'PPR')] + [
        Patch(facecolor='#555555', label='left half: Illumina'), Patch(facecolor='#bbbbbb', label='right half: Aviti'),
        Line2D([], [], color=INK, lw=0.9, label='median of the 23 individuals')]


def draw_support_by_class(axes, HG):
    """HG002: the share of the GIAB true positives, false positives (wrong allele) and false positives (wrong genotype) with k-mer
    support; Illumina on the first axis, Aviti on the second"""
    X = xpos()
    w = 0.27
    for ax, (plat, pl) in zip(axes, PLATS):
        for s in LABEL:
            for k, g in enumerate(('TP', 'FP', 'FP_gt')):
                v = HG[(s, plat)][f'sup_{g}']
                x = X[s] + (k - 1) * w
                ax.bar(x, v, width=w * 0.95, color=CAT[g], edgecolor='white', lw=0.4)
                ax.text(x, v + 1.2, f'{v:.1f}', ha='center', va='bottom', fontsize=5.3, color=INK)
        ax.set_ylim(0, 112)
        ax.set_yticks(range(0, 101, 20))
        ax.grid(axis='y', color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.set_xticks([X[s] for s in LABEL])
        ax.tick_params(axis='x', length=0)
        ax.set_ylabel(f'{pl}\nwith k-mer support (%)')
    axes[1].set_xticklabels([LABEL[s] for s in LABEL], fontsize=6.3)
    return [Patch(color=CAT['TP'], label='GIAB true positives'), Patch(color=CAT['FP'], label='GIAB false positives (allele)'),
            Patch(color=CAT['FP_gt'], label='GIAB false positives (genotype)')]


def draw_region_heatmaps(axes, R):
    """the mean Merfin rate of the cohort by region and set, Illumina and Aviti; returns the image for the colour bar"""
    im = None
    for ax, (plat, pl) in zip(axes, PLATS):
        M = np.array([[np.nanmean(R[(plat, s, st, 'ALL')]) for s in LABEL] for st in ROWS])
        im = ax.imshow(M, cmap='viridis', vmin=40, vmax=100, aspect='auto')
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                ax.text(j, i, f'{M[i, j]:.1f}', ha='center', va='center', fontsize=5.6, color='white' if M[i, j] < 80 else INK)
        ax.set_xticks(range(len(LABEL)))
        ax.set_xticklabels([LABEL[s] for s in LABEL], rotation=90, fontsize=6)
        ax.set_yticks(range(len(ROWS)))
        ax.set_yticklabels([LAB[r] for r in ROWS], fontsize=6.5)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(pl, fontsize=8)
    return im


def draw_counts_by_region(fig, NV, R):
    """per region, the GATK-only (vs PPR) and PPR-only (vs GATK) variants per individual (median; dashed outline) and those with
    k-mer support (filled), with the mean rate above; regions in three groups by size, each with its own axis"""
    sets = ['GATKvsPPR', 'PPRvsGATK']

    def med(plat, s, st, which):
        return float(np.median(NV[(plat, s, st)][which]))
    big = {st: max(med(p, s, st, 0) for p, _ in PLATS for s in sets) for st in REGIONS}
    groups = [g for g in ([st for st in REGIONS if big[st] >= 1e5], [st for st in REGIONS if 1e4 <= big[st] < 1e5],
                          [st for st in REGIONS if big[st] < 1e4]) if g]
    axes = np.atleast_1d(fig.subplots(1, len(groups), gridspec_kw={'width_ratios': [len(g) for g in groups], 'wspace': 0.3}))
    w = 0.21
    for ax, group in zip(axes, groups):
        for i, st in enumerate(group):
            k = 0
            for s in sets:
                for (plat, _), pi in zip(PLATS, (0, 1)):
                    x = i + (k - 1.5) * w
                    tot, val = med(plat, s, st, 0), med(plat, s, st, 1)
                    ax.bar(x, tot, width=w * 0.9, color='none', edgecolor=SHADE[PIPE[s]][pi], lw=0.8, ls='--')
                    ax.bar(x, val, width=w * 0.6, color=SHADE[PIPE[s]][pi], lw=0)
                    ax.text(x, tot * 1.02, f'{np.nanmean(R[(plat, s, st, "ALL")]):.1f}%', ha='center', va='bottom', fontsize=6.6,
                            fontweight='bold', rotation=90, color=SHADE[PIPE[s]][0])
                    k += 1
        ax.set_xticks(range(len(group)))
        ax.set_xticklabels(group, fontsize=7.5)
        ax.tick_params(axis='x', length=0)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.32)
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
            lambda v, _: f'{v / 1000:.0f}k' if v >= 1e4 else f'{v / 1000:.1f}k' if v >= 1000 else f'{v:.0f}'))
        ax.grid(axis='y', color=GRID, lw=0.6)
        ax.set_axisbelow(True)
    axes[0].set_ylabel('Variants per individual (median)')
    return [Patch(color=SHADE['GATK'][0], label='GATK only (vs PPR), Illumina'), Patch(color=SHADE['GATK'][1], label='GATK only, Aviti'),
            Patch(color=SHADE['PPR'][0], label='PPR only (vs GATK), Illumina'), Patch(color=SHADE['PPR'][1], label='PPR only, Aviti'),
            Patch(facecolor='none', edgecolor=INK, ls='--', label='all (dashed); with k-mer support (filled)')]
