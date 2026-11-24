#!/usr/bin/env python3
"""Holm-adjusted p of the genome-level SNV and indel correlations reported in the Results.

coverage_table.py keeps the raw p for these rows. Here they are corrected as the Results report them: the depth difference
against the cross-platform concordance, Holm over the three pipelines; the depth against the concordance of the three pipeline
pairs of one platform, Holm over the three pairs.
Usage: python3 holm_by_variant_type.py [coverage_correlations.tsv]   (writes coverage_correlations_by_type.tsv next to it)
"""
import csv
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'coverage_correlations.tsv')


def holm(ps):
    order, adj, run = sorted(range(len(ps)), key=lambda i: ps[i]), [0.0] * len(ps), 0.0
    for k, i in enumerate(order):
        run = max(run, (len(ps) - k) * ps[i])
        adj[i] = min(1.0, run)
    return adj


fam = defaultdict(list)
with open(SRC, newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['stratum'] != 'genome' or r['type'] not in ('SNV', 'INDEL'):
            continue
        if r['x'] == 'depth difference':
            fam[('depth difference, three pipelines', r['type'])].append(r)
        elif r['family'] == 'within':
            fam[(f"depth, three pairs, {r['platform']}", r['type'])].append(r)

rows = []
for key in sorted(fam):
    for r, a in zip(fam[key], holm([float(r['p']) for r in fam[key]])):
        rows.append([key[0], key[1], r['y'], r['spearman_rho'], r['p'], f'{a:.4g}'])
with open(os.path.join(os.path.dirname(os.path.abspath(SRC)), 'coverage_correlations_by_type.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['family', 'type', 'concordance', 'spearman_rho', 'p', 'p_holm'])
    w.writerows(rows)
for row in rows:
    print('\t'.join(row))
