#!/usr/bin/env python3
"""The data of Supplementary Figure S3 and of the Chromosomes 20-22 sheet of Supplementary File 3: GATK (PASS and ALL), BWA-MEM +
DeepVariant (PASS) and PPR (PASS) on HG002, Illumina and Aviti, scored on chromosomes 20-22 only, in one comparison (DeepVariant holds
these chromosomes out of the training of its models), next to the whole-genome values of Figure 1.
Inputs: <call set>/summary.tsv and BENCH_DONE (happy.sh, region BB), <call set>/chr20_22/CHR20_22_DONE,
<call set>.chr20_22.summary.csv and <call set>.chr20_22.roc.Locations.{SNP,INDEL}.PASS.csv.gz (analyses/benchmark/bench_chr20.sh).
Outputs: chr20_22_benchmark.tsv (precision, recall, F1, TP, FN, FP per platform, pipeline, scope genome or chr20_22, type and filter),
chr20_22_roc.tsv (precision and recall at each QUAL threshold of the PASS calls on chromosomes 20-22) and
chr20_22_data.provenance.txt (md5 of the inputs, of this script and of the tables). It stops if F1 does not match precision and recall
or if the ROC point at the lowest QUAL threshold is not the PASS point of the chromosomes 20-22 summary.
Usage: BENCHMARK_REVISION_DIR=<benchmark folder> python3 chr20_table.py
"""
import csv
import gzip
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BR = os.environ.get('BENCHMARK_REVISION_DIR', '.')
TYPES = ['SNP', 'INDEL']
CALLSETS = [  # platform, pipeline, folder, filters reported (as benchmark_table.py)
    ('Illumina', 'GATK', 'Illumina_GATK_HaplotypeCaller_joint_vqsrbp_snp99.9_indel99.9', ('PASS', 'ALL')),
    ('Illumina', 'BWA-MEM + DeepVariant', 'Illumina_GATK_Deepvariant_markdup', ('PASS',)),
    ('Illumina', 'PPR', 'Illumina_PPR_Deepvariant', ('PASS',)),
    ('Aviti', 'GATK', 'Aviti_GATK_HaplotypeCaller_joint_vqsrbp_snp99.9_indel99.9', ('PASS', 'ALL')),
    ('Aviti', 'BWA-MEM + DeepVariant', 'Aviti_GATK_Deepvariant_markdup', ('PASS',)),
    ('Aviti', 'PPR', 'Aviti_PPR_Deepvariant', ('PASS',)),
]
COLS = ['platform', 'pipeline', 'scope', 'type', 'filter', 'precision', 'recall', 'f1', 'TP', 'FN', 'FP', 'source']
ROC_COLS = ['platform', 'pipeline', 'type', 'qual_threshold', 'precision', 'recall', 'source']


def md5(path):
    with open(path, 'rb') as fh:
        return hashlib.md5(fh.read()).hexdigest()


def fail(msg):
    print('ERROR:', msg)
    print('CHR20_TABLE_FAILED')
    sys.exit(1)


def check_f1(p, rc, f1, what):
    if abs(float(f1) - 2 * float(p) * float(rc) / (float(p) + float(rc))) > 1e-4:
        fail(f'F1 does not match precision and recall: {what}')


rows, roc_rows, prov = [], [], []
num = lambda v: str(int(float(v)))
for plat, pipe, folder, filters in CALLSETS:
    d = os.path.join(BR, folder)
    for marker in ('BENCH_DONE', 'chr20_22/CHR20_22_DONE'):
        if not os.path.exists(os.path.join(d, marker)):
            fail(f'{folder}: no {marker}')
    # whole genome: region BB of summary.tsv (happy.sh)
    sp = os.path.join(d, 'summary.tsv')
    prov.append((sp, md5(sp)))
    recs = list(csv.DictReader(open(sp, newline=''), delimiter='\t'))
    # chromosomes 20-22: hap.py summary.csv (bench_chr20.sh, no strata)
    cp = os.path.join(d, 'chr20_22', f'{folder}.chr20_22.summary.csv')
    prov.append((cp, md5(cp)))
    crecs = list(csv.DictReader(open(cp, newline='')))
    for typ in TYPES:
        for flt in filters:
            hit = [r for r in recs if r['region'] == 'BB' and r['type'] == typ and r['filter'] == flt]
            if len(hit) != 1:
                fail(f'{folder}: {len(hit)} genome rows for {typ} {flt}')
            r = hit[0]
            check_f1(r['precision'], r['recall'], r['F1'], f'{folder} genome {typ} {flt}')
            rows.append([plat, pipe, 'genome', typ, flt, r['precision'], r['recall'], r['F1'], r['TP'], r['FN'], r['FP'],
                         f'{folder}/summary.tsv'])
            hit = [r for r in crecs if r['Type'] == typ and r['Filter'] == flt]
            if len(hit) != 1:
                fail(f'{folder}: {len(hit)} chr20-22 rows for {typ} {flt}')
            r = hit[0]
            check_f1(r['METRIC.Precision'], r['METRIC.Recall'], r['METRIC.F1_Score'], f'{folder} chr20-22 {typ} {flt}')
            rows.append([plat, pipe, 'chr20_22', typ, flt, r['METRIC.Precision'], r['METRIC.Recall'], r['METRIC.F1_Score'],
                         num(r['TRUTH.TP']), num(r['TRUTH.FN']), num(r['QUERY.FP']), f'{folder}/chr20_22/{os.path.basename(cp)}'])
        # precision-recall by QUAL threshold of the PASS calls on chromosomes 20-22 (rows of the whole chromosome: Subtype, Subset and
        # Genotype '*')
        rocf = os.path.join(d, 'chr20_22', f'{folder}.chr20_22.roc.Locations.{typ}.PASS.csv.gz')
        if not os.path.exists(rocf):
            fail(f'{folder}: no {os.path.basename(rocf)}')
        prov.append((rocf, md5(rocf)))
        pts = []
        with gzip.open(rocf, 'rt') as fh:
            for r in csv.DictReader(fh):
                if (r['Type'], r['Subtype'], r['Subset'], r['Genotype'], r['QQ.Field']) != (typ, '*', '*', '*', 'QUAL'):
                    continue
                try:   # thresholds above the highest QUAL have no calls and no precision: left out
                    float(r['METRIC.Precision']), float(r['METRIC.Recall'])
                    pts.append((float(r['QQ']), r['METRIC.Precision'], r['METRIC.Recall']))
                except ValueError:
                    continue
        pts.sort(key=lambda t: t[0])
        if len(pts) < 20:
            fail(f'{rocf}: only {len(pts)} points')
        ps = [x for x in rows if x[:5] == [plat, pipe, 'chr20_22', typ, 'PASS']][0]
        if abs(float(pts[0][1]) - float(ps[5])) > 1e-4 or abs(float(pts[0][2]) - float(ps[6])) > 1e-4:
            fail(f'{folder} {typ}: the ROC point at the lowest QUAL threshold is not the PASS point of the chr20-22 summary')
        roc_rows += [[plat, pipe, typ, f'{q:g}', p, rc, f'{folder}/chr20_22/{os.path.basename(rocf)}'] for q, p, rc in pts]


def write(name, cols, data):
    out = os.path.join(HERE, name)
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(cols)
        w.writerows(data)
    return out


outs = [write('chr20_22_benchmark.tsv', COLS, rows), write('chr20_22_roc.tsv', ROC_COLS, roc_rows)]
with open(os.path.join(HERE, 'chr20_22_data.provenance.txt'), 'w') as fh:
    fh.write(f'chr20_table.py md5 {md5(os.path.abspath(__file__))}\n')
    for p, m in prov:
        fh.write(f'input {m} {p}\n')
    for p in outs:
        fh.write(f'output {md5(p)} {os.path.basename(p)}\n')
for r in rows:
    if r[4] == 'PASS':
        print('\t'.join(r[:8]))
print(f'{len(rows)} rows, {len(roc_rows)} ROC points')
print('CHR20_TABLE_OK')
