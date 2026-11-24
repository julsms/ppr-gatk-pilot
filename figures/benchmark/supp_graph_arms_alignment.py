#!/usr/bin/env python3
"""Supplementary Figure S7: the eight alignment metrics (% of reads) of the five set-ups, Illumina and Aviti together: GAM aligned,
perfect, gapless, properly paired and MAPQ 60 (vg stats -a) and BAM mapped, properly paired and MAPQ 60 (giraffe: after surjection to
chr1-22, X, Y; BWA-MEM: its own BAM, counted on chr1-22, X, Y), drawn with figure2_panels.draw_alignment as Figure 2e. BWA-MEM has no
graph alignment: in the GAM aligned and GAM MAPQ 60 panels its row shows its BAM on all contigs, marked "BAM, all contigs"; in the
other GAM panels its row says "no graph alignment".
Inputs: graph_arms_mapping.tsv. Outputs: figures/SuppFig_graph_arms_alignment.{png,pdf}.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figure2_panels as f2  # noqa: E402
import graph_setups as gs  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

f2.style()
data = f2.load()
fig, axes = plt.subplots(2, 4, figsize=(7.1, 3.6))
f2.draw_alignment(axes.flat, data, metrics=[m for m, _ in f2.METRICS], bwa_in_gam=True)
fig.subplots_adjust(hspace=0.75, wspace=0.12)
fig.legend(handles=f2.setup_handles() + f2.key_handles(('platform',)), loc='upper left', bbox_to_anchor=(0.02, 0.0),
           ncol=1, fontsize=6.3, alignment='left', handlelength=2.6, labelspacing=0.4)
os.makedirs(os.path.join(f2.HERE, 'figures'), exist_ok=True)
for ext, kw in (('png', {'dpi': 600}), ('pdf', {})):
    out = os.path.join(f2.HERE, 'figures', f'SuppFig_graph_arms_alignment.{ext}')
    fig.savefig(out, bbox_inches='tight', facecolor='white', **kw)
    print('written', out)
print('SUPP_GRAPH_ALIGNMENT_DONE')
