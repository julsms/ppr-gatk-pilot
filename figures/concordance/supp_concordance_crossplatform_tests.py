#!/usr/bin/env python3
"""Supplementary Figure S12: concordance across platforms of PPR vs GATK, DV vs GATK and PPR vs DV, mean paired difference by region,
for all variants, SNVs and indels.
Usage: python3 supp_concordance_crossplatform_tests.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure3_panels as f3  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

COMPARISONS = [('H2_PPR_vs_GATK_crossplatform', 'PPR vs GATK'), ('H2_DV_vs_GATK_crossplatform', 'DV vs GATK'),
               ('H2_PPR_vs_DV_crossplatform', 'PPR vs DV')]
f3.style()
D, T = f3.load()
fig, axes = plt.subplots(3, 3, figsize=(7.1, 8.4), sharey=True)
for i, (typ, tl) in enumerate(f3.TYPES):
    for j, (cid, cl) in enumerate(COMPARISONS):
        f3.forest(axes[i, j], T, cid, typ, labels=j == 0)
        axes[i, j].set_title(f'{tl}: {cl}', fontweight='bold')
fig.supxlabel('Difference in genotype concordance across platforms, first minus second pipeline (percentage points)', fontsize=7.5,
              y=0.035)
fig.legend(handles=f3.legend_forest('the first pipeline agrees more across platforms', 'the second pipeline agrees more'),
           loc='upper center', bbox_to_anchor=(0.5, 0.025), ncol=2, fontsize=6.4)
fig.subplots_adjust(wspace=0.42, hspace=0.32, bottom=0.08)
os.makedirs(os.path.join(f3.HERE, 'figures'), exist_ok=True)
for ext, kw in (('png', {'dpi': 600}), ('pdf', {})):
    fig.savefig(os.path.join(f3.HERE, 'figures', f'SuppFig_concordance_crossplatform_tests.{ext}'), bbox_inches='tight',
                facecolor='white', **kw)
print('written figures/SuppFig_concordance_crossplatform_tests')
print('SUPP_CONCORDANCE_S11_DONE')
