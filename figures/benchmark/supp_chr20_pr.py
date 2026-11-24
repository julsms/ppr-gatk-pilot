#!/usr/bin/env python3
"""Supplementary Figure S3: precision and recall of GATK (VQSR 99.9/99.9), BWA-MEM + DeepVariant and PPR on HG002 chromosomes
20-22, which DeepVariant holds out of the training of its models, in one comparison; a SNVs, b indels, Illumina and Aviti
together. One curve per pipeline and platform, precision against recall at each QUAL threshold of the PASS calls (hap.py ROC);
the end of each curve (lowest threshold) is the PASS call set, marked with a circle (Illumina) or a square (Aviti); colours of
Figure 1; Illumina solid lines, Aviti dashed. Drawn at the printed size, 183 x 96 mm.
Inputs: chr20_22_roc.tsv and chr20_22_benchmark.tsv (chr20_table.py); figure1_panels.py (colours, markers, style).
Outputs: figures/SuppFig_chr20_22_PR.pdf, .png and .tiff (600 dpi, TIFF RGB).
"""
import csv
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import figure1_panels as f1  # noqa: E402

LS = {'Illumina': '-', 'Aviti': (0, (3, 2))}
TITLE = {'SNP': 'Chromosomes 20-22, SNVs', 'INDEL': 'Chromosomes 20-22, indels'}
MM = 1 / 25.4

f1.style()
plt.rcParams.update({'pdf.fonttype': 42, 'lines.linewidth': 1.2})
roc = list(csv.DictReader(open(os.path.join(HERE, 'chr20_22_roc.tsv')), delimiter='\t'))
tab = {(r['platform'], r['pipeline'], r['type']): r
       for r in csv.DictReader(open(os.path.join(HERE, 'chr20_22_benchmark.tsv')), delimiter='\t')
       if r['scope'] == 'chr20_22' and r['filter'] == 'PASS'}
if len(tab) != 12:
    raise SystemExit(f'chr20_22_benchmark.tsv: {len(tab)} chr20_22 PASS rows, expected 12')

fig, axes = plt.subplots(1, 2, figsize=(183 * MM, 96 * MM))
fig.subplots_adjust(left=0.085, right=0.985, top=0.9, bottom=0.40, wspace=0.28)
for ax, typ, tag in zip(axes, ('SNP', 'INDEL'), 'ab'):
    ends = [float(tab[(pl, pi, typ)]['recall']) for pl in f1.PLATFORMS for pi in f1.PIPES]
    xlo, xhi = min(ends) - 0.03, max(ends) + 0.006
    shown = []
    for plat in f1.PLATFORMS:
        for pipe in f1.PIPES:
            pts = sorted(((float(r['recall']), float(r['precision'])) for r in roc
                          if (r['platform'], r['pipeline'], r['type']) == (plat, pipe, typ)), key=lambda t: t[0])
            pts = [(x, y) for x, y in pts if x >= xlo]
            shown += [y for _, y in pts]
            c = f1.COLOR[(pipe, plat)]
            ax.plot([x for x, _ in pts], [y for _, y in pts], color=c, ls=LS[plat], lw=1.2, zorder=3)
            e = tab[(plat, pipe, typ)]
            ax.plot(float(e['recall']), float(e['precision']), ls='', marker=f1.MARKER[plat], ms=f1.MSIZE[plat] + 0.5, color=c,
                    mec='white', mew=0.4, zorder=4)
    ylo = min(shown) - (1 - min(shown)) * 0.08
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, 1.0 + (1 - ylo) * 0.02)
    ax.set_xlabel('Recall')
    if typ == 'SNP':
        ax.set_ylabel('Precision')
    ax.grid(color=f1.GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(TITLE[typ], loc='left', fontweight='bold', fontsize=8, pad=8)
    ax.text(-0.12 if typ == 'SNP' else -0.1, 1.06, tag, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')

handles = [Line2D([], [], color=f1.COLOR[(pipe, plat)], ls=LS[plat], lw=1.2, marker=f1.MARKER[plat], ms=f1.MSIZE[plat] + 0.5,
                  mec='white', mew=0.4, label=f'{pipe}, {plat}') for plat in f1.PLATFORMS for pipe in f1.PIPES]
fig.legend(handles=handles, loc='lower center', ncol=2, bbox_to_anchor=(0.5, 0.02), fontsize=7, handlelength=2.6,
           columnspacing=1.6, title='Curves: QUAL thresholds of the PASS calls; marker: the PASS call set', title_fontsize=7)
os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True)
base = os.path.join(HERE, 'figures', 'SuppFig_chr20_22_PR')
fig.savefig(base + '.pdf', facecolor='white')
fig.savefig(base + '.png', dpi=600, facecolor='white')
fig.savefig(base + '.tiff', dpi=600, facecolor='white', pil_kwargs={'compression': 'tiff_lzw'})
try:
    from PIL import Image
    im = Image.open(base + '.tiff')
    if im.mode != 'RGB':
        im.convert('RGB').save(base + '.tiff', compression='tiff_lzw', dpi=(600, 600))
except ImportError:
    pass
print('written', base + '.pdf .png .tiff (183 x 96 mm)')
print('SUPP_CHR20_DONE')
