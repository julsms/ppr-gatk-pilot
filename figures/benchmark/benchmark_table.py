#!/usr/bin/env python3
"""The data of Figure 1 and Supplementary Figures S1 and S2: precision, recall and F1 per region, variant type and filter for GATK (PASS
and ALL), BWA-MEM + DeepVariant (PASS) and PPR (PASS), Illumina and Aviti, copied from <call set>/summary.tsv (CMRG left out).
Outputs: benchmark_data.tsv and benchmark_data.provenance.txt (md5 of the inputs, of this script and of the table).
Usage: BENCHMARK_REVISION_DIR=<benchmark folder> python3 benchmark_table.py
"""
import csv
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BR = os.environ.get('BENCHMARK_REVISION_DIR', '.')
CALLSETS = [  # platform, pipeline, folder, filters reported
    ('Illumina', 'GATK', 'Illumina_GATK_HaplotypeCaller_joint_vqsrbp_snp99.9_indel99.9', ('PASS', 'ALL')),
    ('Illumina', 'BWA-MEM + DeepVariant', 'Illumina_GATK_Deepvariant_markdup', ('PASS',)),
    ('Illumina', 'PPR', 'Illumina_PPR_Deepvariant', ('PASS',)),
    ('Aviti', 'GATK', 'Aviti_GATK_HaplotypeCaller_joint_vqsrbp_snp99.9_indel99.9', ('PASS', 'ALL')),
    ('Aviti', 'BWA-MEM + DeepVariant', 'Aviti_GATK_Deepvariant_markdup', ('PASS',)),
    ('Aviti', 'PPR', 'Aviti_PPR_Deepvariant', ('PASS',)),
]
REGIONS = ['BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']   # CMRG left out
TYPES = ['SNP', 'INDEL']
COLS = ['platform', 'pipeline', 'region', 'type', 'filter', 'precision', 'recall', 'f1', 'TP', 'FN', 'FP', 'truth_total',
        'query_total', 'source']


def md5(path):
    with open(path, 'rb') as fh:
        return hashlib.md5(fh.read()).hexdigest()


def fail(msg):
    print('ERROR:', msg)
    print('BENCHMARK_TABLE_FAILED')
    sys.exit(1)


rows, prov = [], []
for plat, pipe, folder, filters in CALLSETS:
    d = os.path.join(BR, folder)
    summ = os.path.join(d, 'summary.tsv')
    if not os.path.exists(os.path.join(d, 'BENCH_DONE')) or not os.path.exists(summ):
        fail(f'{folder}: no BENCH_DONE or summary.tsv')
    prov.append((summ, md5(summ)))
    with open(summ, newline='') as fh:
        recs = list(csv.DictReader(fh, delimiter='\t'))
    for reg in REGIONS:
        for typ in TYPES:
            for flt in filters:
                hit = [r for r in recs if r['region'] == reg and r['type'] == typ and r['filter'] == flt]
                if len(hit) != 1:
                    fail(f'{folder}: {len(hit)} rows for {reg} {typ} {flt}')
                r = hit[0]
                p, rc, f1 = float(r['precision']), float(r['recall']), float(r['F1'])
                if abs(f1 - 2 * p * rc / (p + rc)) > 1e-4:
                    fail(f'{folder}: F1 does not match precision and recall for {reg} {typ} {flt}')
                rows.append([plat, pipe, reg, typ, flt, r['precision'], r['recall'], r['F1'], r['TP'], r['FN'], r['FP'],
                             r['truth_total'], r['query_total'], f'{folder}/summary.tsv'])

expected = 2 * (2 + 1 + 1) * len(REGIONS) * len(TYPES)
if len(rows) != expected:
    fail(f'{len(rows)} rows, expected {expected}')
out = os.path.join(HERE, 'benchmark_data.tsv')
with open(out, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(COLS)
    w.writerows(rows)
with open(os.path.join(HERE, 'benchmark_data.provenance.txt'), 'w') as fh:
    for path, m in prov:
        fh.write(f'input\t{m}\t{path}\n')
    fh.write(f'script\t{md5(os.path.abspath(__file__))}\t{os.path.abspath(__file__)}\n')
    fh.write(f'table\t{md5(out)}\t{out}\n')
print(f'{len(rows)} rows written to {out}')
print('BENCHMARK_TABLE_OK')
