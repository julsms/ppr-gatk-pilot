#!/usr/bin/env python3
"""supp_uniquevariants_consequence_classes.py - Supplementary Figure S16: the predicted consequences of the pipeline-unique variants
(VEP, most severe consequence). a the share of the distinct alleles of the cohort in five groups, genome, with the share of coding or
splice alleles at the end of each bar; b the coding and splice consequences. Input: uv_consequences.tsv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure4_panels as f4  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

GROUPS = [('intergenic', '#DDCC77', '#1f2933'), ('intron', '#88CCEE', '#1f2933'),
          ('upstream or downstream of a gene', '#44AA99', '#1f2933'),
          ('UTR, non-coding exon, regulatory or other', '#332288', '#ffffff'), ('coding or splice', '#CC3366', '#ffffff')]
SETS = [('GATK', 'PPR'), ('PPR', 'GATK'), ('DV', 'PPR'), ('PPR', 'DV'), ('PPR', 'BOTH')]
f4.style()
n = {(r['set'], r['platform'], r['group']): (int(r['alleles']), int(r['alleles_in_region'])) for r in f4.read('uv_consequences.tsv')
     if r['region'] == 'genome' and r['kind'] == 'composition'}
_, _, Q = f4.load()
fig = plt.figure(figsize=(7.1, 6.6))
gs = fig.add_gridspec(4, 1, height_ratios=[3.0, 0.3, 2.0, 0.4], hspace=0.45)
# a: the five groups
ax = fig.add_subplot(gs[0])
labels, ypos, y = [], [], 0.0
for a, b in SETS:
    for pi, (plat, pl) in enumerate(f4.PLATFORMS):
        left = 0.0
        total = n[(f'{a}vs{b}', plat, GROUPS[0][0])][1]
        for g, col, ink in GROUPS:
            v, tot = n[(f'{a}vs{b}', plat, g)]
            w = 100 * v / tot
            ax.barh(y, w, left=left, height=0.8, color=col, edgecolor='white', lw=0.4)
            if w >= 6:
                ax.text(left + w / 2, y, f'{w:.0f}', ha='center', va='center', fontsize=6, color=ink)
            left += w
        if abs(left - 100) > 1e-6:
            raise SystemExit(f'{a}vs{b} {plat}: the groups add up to {left:.6f} %')
        coding = 100 * n[(f'{a}vs{b}', plat, 'coding or splice')][0] / total
        ax.text(101, y, f'{coding:.1f}', ha='left', va='center', fontsize=5.8, color='#CC3366', fontweight='bold')
        labels.append(f'{f4.set_label(a, b)}, {pl}')
        ypos.append(y)
        y -= 1.0
    y -= 0.5
ax.set_yticks(ypos)
ax.set_yticklabels(labels, fontsize=6.3)
ax.tick_params(axis='y', length=0)
ax.set_xlim(0, 100)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel('Share of the pipeline-unique variants of the cohort, genome (%)')
ax.text(101, max(ypos) + 0.9, 'coding or\nsplice (%)', ha='left', va='bottom', fontsize=5.8, color='#CC3366')
lax = fig.add_subplot(gs[1])
lax.axis('off')
lax.legend(handles=[Patch(color=col, label=g) for g, col, _ in GROUPS], loc='center', ncol=3, fontsize=6.3)
# b: coding and splice consequences
bx = fig.add_subplot(gs[2])
f4.draw_consequences(bx, Q)
lbx = fig.add_subplot(gs[3])
lbx.axis('off')
lbx.legend(handles=f4.legend_consequences(), loc='center', ncol=4, fontsize=5.8, handlelength=1.6, columnspacing=1.0)
f4.panel_tags(fig, [(ax, 'a', 'Consequence groups'), (bx, 'b', 'Coding and splice consequences')])
f4.supp_title(fig, 'Predicted consequences of the pipeline-unique variants')
f4.save(fig, 'SuppFig_uniquevariants_consequence_classes')
print('SUPP_UV_CONSEQUENCES_DONE')
