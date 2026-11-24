#!/usr/bin/env python3
"""supp_kmer_spectra.py - supplementary figure: k-mer spectra (k = 21, FastK histograms of the depth-matched reads) of every individual
on both platforms, from GENOMESCOPE_DIR (default genomescope): sample_pairs.csv and genomescope_results_{illumina,aviti}/. The last,
cumulative bin is dropped. Output: figures/SuppFig_kmer_spectra.{pdf,png}.
"""
import hashlib
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FuncFormatter, MaxNLocator  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
GS = os.environ.get('GENOMESCOPE_DIR', 'genomescope')
DPI = int(os.environ.get('DPI', '600'))
COL = {'Illumina': '#2E86AB', 'Aviti': '#A23B72'}
INK, GRID = '#1f2933', '#e3e7eb'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 6.5, 'axes.labelsize': 7, 'xtick.labelsize': 5.6,
                     'ytick.labelsize': 5.6, 'axes.edgecolor': '#9aa5b1', 'axes.labelcolor': INK, 'xtick.color': INK,
                     'ytick.color': INK, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
                     'axes.unicode_minus': False})


def fail(msg):
    print('ERROR:', msg)
    print('SUPP_KMER_SPECTRA_FAILED')
    sys.exit(1)


def md5(path):
    with open(path, 'rb') as fh:
        return hashlib.md5(fh.read()).hexdigest()


def pairs():
    """(Illumina number, illumina id, aviti id) from the pairs file, in increasing Illumina number"""
    out = []
    with open(os.path.join(GS, 'sample_pairs.csv')) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            ill, av = [x.strip() for x in line.split(',')]
            out.append((int(re.search(r'\d+', ill).group()), ill, av))
    out.sort()
    by = {n: int(re.search(r'\d+', av).group()) for n, _, av in out}
    if by.get(12) != 21 or by.get(21) != 12:
        fail('the pairs file is not the sample_pairs.csv of analyses/genomescope')
    if len(out) != 23:
        fail(f'{len(out)} pairs, expected 23')
    return out


def histo(platform_dir, sample_id):
    """the k = 21 histogram of one sample: its folder starts with the sample id; the last bin is cumulative and is dropped"""
    root = os.path.join(GS, platform_dir)
    dirs = [d for d in os.listdir(root) if d == sample_id or d.startswith(sample_id + '_')]
    if len(dirs) != 1:
        fail(f'{platform_dir}: {len(dirs)} folders for {sample_id}')
    files = sorted(f for f in os.listdir(os.path.join(root, dirs[0])) if f.endswith('_k21.histo'))
    if len(files) != 1:
        fail(f'{platform_dir}/{dirs[0]}: {len(files)} k21 histograms')
    path = os.path.join(root, dirs[0], files[0])
    rows = []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                rows.append((int(parts[0]), int(parts[1])))
    if len(rows) < 10:
        fail(f'{path}: {len(rows)} rows')
    a = np.array(rows[:-1], dtype=float)
    return a[:, 0], a[:, 1], path


def peak(cov, freq, min_cov=4):
    """coverage and frequency of the homozygous peak: the highest bin above the error k-mers (coverage above min_cov)"""
    keep = cov > min_cov
    i = int(np.argmax(freq[keep]))
    return cov[keep][i], freq[keep][i]


P = pairs()
data, paths = [], []
for k, (n, ill, av) in enumerate(P, start=1):
    ci, fi, pi = histo('genomescope_results_illumina', ill)
    ca, fa, pa = histo('genomescope_results_aviti', av)
    data.append((k, (ci, fi), (ca, fa)))
    paths += [pi, pa]
peaks = [peak(*h) for _, hi, ha in data for h in (hi, ha)]
x_max = int(max(40, min(150, np.median([c for c, _ in peaks]) * 4)))

W, ncol = 7.087, 4
nrow = -(-len(data) // ncol)
fig, axes = plt.subplots(nrow, ncol, figsize=(W, 1.32 * nrow + 0.45), sharex=True)
for ax, (k, (ci, fi), (ca, fa)) in zip(axes.flat, data):
    cap = 1.5 * max(peak(ci, fi)[1], peak(ca, fa)[1])
    for (cov, freq), plat in (((ci, fi), 'Illumina'), ((ca, fa), 'Aviti')):
        sel = (cov >= 1) & (cov <= x_max)
        ax.plot(cov[sel], freq[sel], color=COL[plat], lw=1.0)
    ax.set_xlim(0, x_max)
    ax.set_ylim(0, cap)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3, integer=False))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v / 1e6:.0f}'))
    ax.grid(color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.set_title(f'Individual {k}', fontsize=6.5, pad=3)
    ax.tick_params(axis='x', labelbottom=True)
for ax in list(axes.flat)[len(data):]:
    ax.set_visible(False)
fig.subplots_adjust(left=0.08, right=0.99, top=1 - 0.30 / fig.get_figheight(), bottom=0.62 / fig.get_figheight(), hspace=0.55,
                    wspace=0.32)
fig.supylabel('K-mers (millions)', fontsize=7, x=0.01)
fig.supxlabel('K-mer coverage', fontsize=7, y=0.30 / fig.get_figheight())
fig.legend(handles=[Line2D([], [], color=COL[p], lw=1.6, label=p) for p in ('Illumina', 'Aviti')], loc='lower center', ncol=2,
           fontsize=7, bbox_to_anchor=(0.5, 0.0), borderaxespad=0.2)
fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
fig.text(bb.x0 / fig.get_figwidth(), (bb.y1 + 0.08) / fig.get_figheight(), 'K-mer spectra of the paired samples', ha='left',
         va='bottom', fontsize=8.5, fontweight='bold', color=INK)
bad = [t.get_text() for t in fig.findobj(plt.Text) if any(ord(ch) in (0x2012, 0x2013, 0x2014, 0x2015, 0x2212) for ch in t.get_text())]
if bad:
    fail(f'long dash in the figure text: {bad[:3]}')
os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True)
for ext, kw in (('pdf', {}), ('png', {'dpi': DPI})):
    fig.savefig(os.path.join(HERE, 'figures', f'SuppFig_kmer_spectra.{ext}'), bbox_inches='tight', pad_inches=0.04, facecolor='white',
                **kw)
plt.close(fig)
finger = hashlib.md5(''.join(md5(p) for p in paths).encode()).hexdigest()
with open(os.path.join(HERE, 'supp_kmer_spectra.provenance.txt'), 'w') as fh:
    fh.write(f'script\t{md5(os.path.abspath(__file__))}\t{os.path.abspath(__file__)}\n')
    fh.write(f'pairs\t{md5(os.path.join(GS, "sample_pairs.csv"))}\t{os.path.join(GS, "sample_pairs.csv")}\n')
    fh.write(f'histograms\t{finger}\t{len(paths)} k21 histograms under {GS} (md5 of their md5 sums, in the order of the panels)\n')
    fh.write(f'matplotlib\t{matplotlib.__version__}\n')
print(f'written figures/SuppFig_kmer_spectra.pdf and .png ({DPI} dpi); {len(data)} individuals, x up to {x_max}')
print('SUPP_KMER_SPECTRA_DONE')
