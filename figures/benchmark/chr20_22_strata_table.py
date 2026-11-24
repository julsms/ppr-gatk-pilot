#!/usr/bin/env python3
"""The HG002 benchmark of Figures 1 and 2 by region on chromosomes 20-22 only (DeepVariant holds these chromosomes out of the training
of its models): the tables of Figure 1 (benchmark_data.tsv) and Figure 2 (graph_arms_data.tsv, graph_arms_roc.tsv) in the same layout,
made from the chromosome 20-22 hap.py runs with strata (analyses/benchmark/bench_chr20.sh with STRATA=1). The sheet "Chromosomes 20-22
by region" of Supplementary File 3 holds these values. graph_arms_mapping.tsv (alignment metrics of the whole genome) is copied
unchanged, so that the figure scripts of this folder can draw the chromosome 20-22 version from the output folder.
Inputs: <call set>/chr20_22_strata/{CHR20_22_STRATA_DONE, <call set>.chr20_22_strata.extended.csv,
<call set>.chr20_22_strata.roc.Locations.{SNP,INDEL}.csv.gz} of the 12 call sets (the six of Figure 1 and the three other graph
set-ups on both platforms); <call set>/chr20_22/<call set>.chr20_22.summary.csv of the six Figure 1 call sets (check);
graph_arms_mapping.tsv of this folder (graph_arms_table.py).
Outputs: chr20_22_strata/{benchmark_data.tsv, graph_arms_data.tsv, graph_arms_roc.tsv, graph_arms_mapping.tsv, provenance.txt}; last
line CHR20_22_STRATA_TABLE_OK. It stops if F1 does not match precision and recall, if a region, type or filter is missing, if the BB
rows differ from the run without strata, or if the ROC point at the lowest QUAL threshold is not the PASS point. ROC points with no
calls left (hap.py leaves precision empty) are left out and counted.
Usage: BENCHMARK_REVISION_DIR=<benchmark folder> python3 chr20_22_strata_table.py
"""
import csv
import gzip
import hashlib
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import graph_setups as gs  # noqa: E402

BR = os.environ.get('BENCHMARK_REVISION_DIR', '.')
TAG = 'chr20_22_strata'
OUT = os.path.join(HERE, TAG)
TYPES = ['SNP', 'INDEL']
REGIONS = gs.REGIONS
FIG1 = [  # as benchmark_table.py
    ('Illumina', 'GATK', 'Illumina_GATK_HaplotypeCaller_joint_vqsrbp_snp99.9_indel99.9', ('PASS', 'ALL')),
    ('Illumina', 'BWA-MEM + DeepVariant', 'Illumina_GATK_Deepvariant_markdup', ('PASS',)),
    ('Illumina', 'PPR', 'Illumina_PPR_Deepvariant', ('PASS',)),
    ('Aviti', 'GATK', 'Aviti_GATK_HaplotypeCaller_joint_vqsrbp_snp99.9_indel99.9', ('PASS', 'ALL')),
    ('Aviti', 'BWA-MEM + DeepVariant', 'Aviti_GATK_Deepvariant_markdup', ('PASS',)),
    ('Aviti', 'PPR', 'Aviti_PPR_Deepvariant', ('PASS',)),
]
num = lambda v: str(int(float(v)))


def md5(path):
    with open(path, 'rb') as fh:
        return hashlib.md5(fh.read()).hexdigest()


def fail(msg):
    print('ERROR:', msg)
    print('CHR20_22_STRATA_TABLE_FAILED')
    sys.exit(1)


prov = []


def extended(folder):
    d = os.path.join(BR, folder, TAG)
    if not os.path.exists(os.path.join(d, 'CHR20_22_STRATA_DONE')):
        fail(f'{folder}: no {TAG}/CHR20_22_STRATA_DONE')
    p = os.path.join(d, f'{folder}.{TAG}.extended.csv')
    prov.append((p, md5(p)))
    ext = {}
    for r in csv.DictReader(open(p, newline='')):
        if r['Type'] in TYPES and r['Subtype'] == '*' and r['Genotype'] == '*' and r['Filter'] in ('ALL', 'PASS'):
            k = ('BB' if r['Subset'] == '*' else r['Subset'], r['Type'], r['Filter'])
            if k in ext:
                fail(f'{p}: two rows for {k}')
            ext[k] = r
    for reg in REGIONS:
        for typ in TYPES:
            for flt in ('ALL', 'PASS'):
                if (reg, typ, flt) not in ext:
                    fail(f'{p}: no row {reg} {typ} {flt}')
                r = ext[(reg, typ, flt)]
                pr, rc, f1 = float(r['METRIC.Precision']), float(r['METRIC.Recall']), float(r['METRIC.F1_Score'])
                if abs(f1 - 2 * pr * rc / (pr + rc)) > 1e-4:
                    fail(f'{p}: F1 does not match precision and recall for {reg} {typ} {flt}')
    return ext, f'{folder}/{TAG}/{os.path.basename(p)}'


os.makedirs(OUT, exist_ok=True)
# Figure 1
rows1 = []
for plat, pipe, folder, filters in FIG1:
    ext, src = extended(folder)
    sp = os.path.join(BR, folder, 'chr20_22', f'{folder}.chr20_22.summary.csv')
    prov.append((sp, md5(sp)))
    for r in csv.DictReader(open(sp, newline='')):
        if r['Type'] in TYPES and r['Filter'] in ('ALL', 'PASS'):
            e = ext[('BB', r['Type'], r['Filter'])]
            for a, b in (('TRUTH.TP', 'TRUTH.TP'), ('TRUTH.FN', 'TRUTH.FN'), ('QUERY.FP', 'QUERY.FP'),
                         ('METRIC.F1_Score', 'METRIC.F1_Score')):
                if abs(float(r[a]) - float(e[b])) > 1e-9:
                    fail(f'{folder} BB {r["Type"]} {r["Filter"]}: {a} {r[a]} without strata, {e[b]} with strata')
    for reg in REGIONS:
        for typ in TYPES:
            for flt in filters:
                r = ext[(reg, typ, flt)]
                rows1.append([plat, pipe, reg, typ, flt, r['METRIC.Precision'], r['METRIC.Recall'], r['METRIC.F1_Score'],
                              num(r['TRUTH.TP']), num(r['TRUTH.FN']), num(r['QUERY.FP']), num(r['TRUTH.TOTAL']),
                              num(r['QUERY.TOTAL']), src])
print('check: the BB rows of the six Figure 1 call sets equal the chromosome 20-22 run without strata')
# Figure 2
rows2, roc = [], []
dropped = 0
for plat in gs.PLATFORMS:
    for key in gs.SETUPS:
        folder = f'{plat}_{gs.CALL_FOLDER[key]}'
        ext, src = extended(folder)
        for reg in REGIONS:
            for typ in TYPES:
                r = ext[(reg, typ, 'PASS')]
                rows2.append([plat, key, reg, typ, r['METRIC.Precision'], r['METRIC.Recall'], r['METRIC.F1_Score'],
                              num(r['TRUTH.TP']), num(r['TRUTH.FN']), num(r['QUERY.FP']), src])
        for typ in TYPES:
            rocf = os.path.join(BR, folder, TAG, f'{folder}.{TAG}.roc.Locations.{typ}.csv.gz')
            prov.append((rocf, md5(rocf)))
            pts = []
            with gzip.open(rocf, 'rt') as fh:
                for r in csv.DictReader(fh):
                    if (r['Type'], r['Subtype'], r['Subset'], r['Filter'], r['Genotype'], r['QQ.Field']) != (typ, '*', '*', 'ALL', '*', 'QUAL'):
                        fail(f'{rocf}: unexpected row {r["Type"]} {r["Subtype"]} {r["Subset"]} {r["Filter"]} {r["Genotype"]}')
                    # at the top QUAL thresholds no calls are left on chromosomes 20-22 and hap.py leaves precision empty (0/0)
                    if r['METRIC.Precision'] == '':
                        dropped += 1
                        continue
                    try:
                        pts.append((float(r['QQ']), r['METRIC.Precision'], r['METRIC.Recall']))
                    except ValueError:
                        continue
            pts.sort(key=lambda t: t[0])
            if len(pts) < 20:
                fail(f'{rocf}: only {len(pts)} points')
            ps = ext[('BB', typ, 'PASS')]
            if abs(float(pts[0][1]) - float(ps['METRIC.Precision'])) > 1e-4 or abs(float(pts[0][2]) - float(ps['METRIC.Recall'])) > 1e-4:
                fail(f'{folder} {typ}: the ROC point at the lowest QUAL threshold is not the PASS point')
            roc += [[plat, key, typ, f'{q:g}', p, rc, f'{folder}/{TAG}/{os.path.basename(rocf)}'] for q, p, rc in pts]


def write(name, cols, rows):
    out = os.path.join(OUT, name)
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(cols)
        w.writerows(rows)
    print(f'{len(rows)} rows written to {out}')
    return out


if len(rows1) != 2 * 4 * len(REGIONS) * len(TYPES) or len(rows2) != 2 * len(gs.SETUPS) * len(REGIONS) * len(TYPES):
    fail(f'{len(rows1)} Figure 1 rows, {len(rows2)} Figure 2 rows')
outs = [write('benchmark_data.tsv', ['platform', 'pipeline', 'region', 'type', 'filter', 'precision', 'recall', 'f1', 'TP', 'FN', 'FP',
                                     'truth_total', 'query_total', 'source'], rows1),
        write('graph_arms_data.tsv', ['platform', 'setup', 'region', 'type', 'precision', 'recall', 'f1', 'TP', 'FN', 'FP', 'source'],
              rows2),
        write('graph_arms_roc.tsv', ['platform', 'setup', 'type', 'qual_threshold', 'precision', 'recall', 'source'], roc)]
shutil.copyfile(os.path.join(HERE, 'graph_arms_mapping.tsv'), os.path.join(OUT, 'graph_arms_mapping.tsv'))
outs.append(os.path.join(OUT, 'graph_arms_mapping.tsv'))
with open(os.path.join(OUT, 'provenance.txt'), 'w') as fh:
    for path, m in prov:
        fh.write(f'input\t{m}\t{path}\n')
    fh.write(f'script\t{md5(os.path.abspath(__file__))}\t{os.path.abspath(__file__)}\n')
    for o in outs:
        fh.write(f'table\t{md5(o)}\t{o}\n')
print(f'ROC points without calls (precision empty) left out: {dropped}')
print('CHR20_22_STRATA_TABLE_OK')
