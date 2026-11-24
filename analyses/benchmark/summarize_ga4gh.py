#!/usr/bin/env python3
"""One table per call set from the hap.py outputs of happy.sh: BB (the whole v5.0q benchmark bed), the strata and CMRG_v1.00_truth,
for SNP and INDEL, ALL and PASS. Columns: region, type, filter, recall, precision, F1, TP, FN, FP, truth_total, query_total, query_unk.
Usage: python3 summarize_ga4gh.py <folder> <name> <strata_ga4gh.with_chr.tsv> > <folder>/summary.tsv
"""
import csv
import os
import sys

D, NAME = sys.argv[1], sys.argv[2]
STRATA = [l.split('\t')[0] for l in open(sys.argv[3]) if l.strip()]
COLS = ['region', 'type', 'filter', 'recall', 'precision', 'F1', 'TP', 'FN', 'FP', 'truth_total', 'query_total', 'query_unk']
num = lambda v: str(int(float(v))) if v not in ('', None) else ''


def row(region, r):
    return [region, r['Type'], r['Filter'], r['METRIC.Recall'], r['METRIC.Precision'], r['METRIC.F1_Score'], num(r['TRUTH.TP']),
            num(r['TRUTH.FN']), num(r['QUERY.FP']), num(r['TRUTH.TOTAL']), num(r['QUERY.TOTAL']), num(r['QUERY.UNK'])]


ext = {}
for r in csv.DictReader(open(os.path.join(D, 'happy', f'{NAME}.extended.csv'))):
    if r['Type'] in ('SNP', 'INDEL') and r['Subtype'] == '*' and r['Genotype'] == '*' and r['Filter'] in ('ALL', 'PASS'):
        ext[('BB' if r['Subset'] == '*' else r['Subset'], r['Type'], r['Filter'])] = r
print('\t'.join(COLS))
for reg in ['BB'] + STRATA:
    for t in ('SNP', 'INDEL'):
        for f in ('ALL', 'PASS'):
            if (reg, t, f) not in ext:
                sys.exit(f'missing in extended.csv: {reg} {t} {f}')
            print('\t'.join(row(reg, ext[(reg, t, f)])))
for r in csv.DictReader(open(os.path.join(D, 'cmrg', f'{NAME}.summary.csv'))):
    if r['Type'] in ('SNP', 'INDEL') and r['Filter'] in ('ALL', 'PASS'):
        print('\t'.join(row('CMRG_v1.00_truth', r)))
