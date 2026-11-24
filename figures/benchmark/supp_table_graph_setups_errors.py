#!/usr/bin/env python3
"""Supplementary Table S1: errors of the five set-ups of the graph decomposition on HG002, PASS calls in the GIAB v5.0q benchmark
regions (BB): false negatives (FN) and false positives (FP) of SNVs and indels, total errors (their sum) and the change of the total
against BWA-MEM in percent, per platform. These are the changes against BWA-MEM in the Results; the change of PPR against HPRC default
(the abstract) is printed too.
Inputs: graph_arms_data.tsv (graph_arms_table.py, same folder). Outputs: graph_setups_errors.tsv; stdout; last line
SUPP_TABLE_GRAPH_SETUPS_DONE.
Usage: python3 supp_table_graph_setups_errors.py
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import graph_setups as gs  # noqa: E402

with open(os.path.join(HERE, 'graph_arms_data.tsv'), newline='') as fh:
    bb = {(r['platform'], r['setup'], r['type']): r for r in csv.DictReader(fh, delimiter='\t') if r['region'] == 'BB'}
if len(bb) != len(gs.PLATFORMS) * len(gs.SETUPS) * 2:
    sys.exit(f'{len(bb)} BB rows in graph_arms_data.tsv, expected {len(gs.PLATFORMS) * len(gs.SETUPS) * 2}')
rows, total = [], {}
for p in gs.PLATFORMS:
    for s in gs.SETUPS:
        v = [int(bb[p, s, t][k]) for t in ('SNP', 'INDEL') for k in ('FN', 'FP')]
        total[p, s] = sum(v)
        change = 100 * (total[p, s] - total[p, 'bwa_grch38']) / total[p, 'bwa_grch38']
        rows.append([p, gs.SHORT[s], *v, total[p, s], '-' if s == 'bwa_grch38' else f'{change:+.1f}'])
        print(f'{p} {gs.SHORT[s]}: SNV FN {v[0]:,} FP {v[1]:,}, indel FN {v[2]:,} FP {v[3]:,}, total {total[p, s]:,}, '
              f'vs BWA-MEM {rows[-1][-1]}')
    ppr_vs_default = 100 * (total[p, 'ppr'] - total[p, 'giraffe_hprc_default']) / total[p, 'giraffe_hprc_default']
    print(f'{p}: PPR vs HPRC default {ppr_vs_default:+.1f} %')
with open(os.path.join(HERE, 'graph_setups_errors.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['platform', 'setup', 'snv_FN', 'snv_FP', 'indel_FN', 'indel_FP', 'total_errors', 'change_vs_bwa_mem_pct'])
    w.writerows(rows)
print('SUPP_TABLE_GRAPH_SETUPS_DONE')
