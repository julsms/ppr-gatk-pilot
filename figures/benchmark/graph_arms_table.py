#!/usr/bin/env python3
"""The data of Figure 2 and Supplementary Figures S5-S7: PASS precision, recall, F1, TP, FN and FP of the five set-ups per platform, region
and type (graph_arms_data.tsv), precision and recall by QUAL threshold from the hap.py ROC files (graph_arms_roc.tsv), the alignment
metrics of the four giraffe set-ups (graph_arms_mapping.tsv, from vg stats and mapping_metrics.tsv) and the BAM metrics of BWA-MEM
(bam_metrics.sh of analyses/benchmark/decomposition on the BWA-MEM BAM: bam_* counted on chr1-22, X, Y in bam_chr1_22XY/, bam_all_* on
all contigs in bam/), and their provenance (md5).
It stops if F1 does not match precision and recall, the lowest-QUAL ROC point is not the PASS point of summary.tsv, or the BWA-MEM BAM
does not have the number of primary records of the giraffe BAMs (the same reads).
Usage: BENCHMARK_REVISION_DIR=<benchmark folder> python3 graph_arms_table.py
"""
import csv
import gzip
import hashlib
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import graph_setups as gs  # noqa: E402

BR = os.environ.get('BENCHMARK_REVISION_DIR', '.')
TYPES = ['SNP', 'INDEL']
COLS = ['platform', 'setup', 'region', 'type', 'precision', 'recall', 'f1', 'TP', 'FN', 'FP', 'source']
ROC_COLS = ['platform', 'setup', 'type', 'qual_threshold', 'precision', 'recall', 'source']
MAP_COLS = ['platform', 'setup', 'metric', 'percent', 'source']
# vg stats -a line -> metric; the percentages are of "Total alignments" (one alignment per read: no secondary alignments)
VG = [('gam_aligned', r'Total aligned: (\d+)', 'graph_aligned_pct'), ('gam_perfect', r'Total perfect: (\d+)', 'graph_perfect_pct'),
      ('gam_gapless', r'Total gapless \(softclips allowed\): (\d+)', 'graph_gapless_pct'),
      ('gam_properly_paired', r'Total properly paired: (\d+)', 'graph_properly_paired_pct'),
      ('gam_mapq60', r'Mapping quality: .*max 60 \((\d+) reads\)', 'graph_mapq60_pct')]
BAM = [('bam_mapped', 'bam_mapped_pct'), ('bam_properly_paired', 'bam_properly_paired_pct'), ('bam_mapq60', 'bam_mapq60_pct')]


def md5(path):
    with open(path, 'rb') as fh:
        return hashlib.md5(fh.read()).hexdigest()


def fail(msg):
    print('ERROR:', msg)
    print('GRAPH_ARMS_TABLE_FAILED')
    sys.exit(1)


def write(name, cols, rows):
    out = os.path.join(HERE, name)
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(cols)
        w.writerows(rows)
    print(f'{len(rows)} rows written to {out}')
    return out


rows, roc_rows, map_rows, prov = [], [], [], []
for plat in gs.PLATFORMS:
    for key in gs.SETUPS:
        folder = f'{plat}_{gs.CALL_FOLDER[key]}'
        summ = os.path.join(BR, folder, 'summary.tsv')
        if not os.path.exists(os.path.join(BR, folder, 'BENCH_DONE')) or not os.path.exists(summ):
            fail(f'{folder}: no BENCH_DONE or summary.tsv')
        prov.append((summ, md5(summ)))
        recs = list(csv.DictReader(open(summ, newline=''), delimiter='\t'))
        for reg in gs.REGIONS:
            for typ in TYPES:
                hit = [r for r in recs if r['region'] == reg and r['type'] == typ and r['filter'] == 'PASS']
                if len(hit) != 1:
                    fail(f'{folder}: {len(hit)} rows for {reg} {typ} PASS')
                r = hit[0]
                p, rc, f1 = float(r['precision']), float(r['recall']), float(r['F1'])
                if abs(f1 - 2 * p * rc / (p + rc)) > 1e-4:
                    fail(f'{folder}: F1 does not match precision and recall for {reg} {typ}')
                rows.append([plat, key, reg, typ, r['precision'], r['recall'], r['F1'], r['TP'], r['FN'], r['FP'],
                             f'{folder}/summary.tsv'])
        # precision-recall by QUAL threshold: the rows of the whole benchmark (Subtype, Subset and Genotype '*'), all records
        for typ in TYPES:
            rocf = os.path.join(BR, folder, 'happy', f'{folder}.roc.Locations.{typ}.csv.gz')
            if not os.path.exists(rocf):
                fail(f'{folder}: no {os.path.basename(rocf)}')
            prov.append((rocf, md5(rocf)))
            pts = []
            with gzip.open(rocf, 'rt') as fh:
                for r in csv.DictReader(fh):
                    if (r['Type'], r['Subtype'], r['Subset'], r['Filter'], r['Genotype'], r['QQ.Field']) != (typ, '*', '*', 'ALL', '*', 'QUAL'):
                        fail(f'{rocf}: unexpected row {r["Type"]} {r["Subtype"]} {r["Subset"]} {r["Filter"]} {r["Genotype"]}')
                    try:
                        pts.append((float(r['QQ']), r['METRIC.Precision'], r['METRIC.Recall']))
                    except ValueError:
                        continue
            pts.sort(key=lambda t: t[0])
            if len(pts) < 20:
                fail(f'{rocf}: only {len(pts)} points')
            ps = [r for r in recs if r['region'] == 'BB' and r['type'] == typ and r['filter'] == 'PASS'][0]
            if abs(float(pts[0][1]) - float(ps['precision'])) > 1e-4 or abs(float(pts[0][2]) - float(ps['recall'])) > 1e-4:
                fail(f'{folder} {typ}: the ROC point at the lowest QUAL threshold is not the PASS point of summary.tsv')
            roc_rows += [[plat, key, typ, f'{q:g}', p, rc, f'{folder}/happy/{os.path.basename(rocf)}'] for q, p, rc in pts]
    # alignment metrics of the giraffe set-ups: GAM from vg stats -a, BAM (after surjection) from mapping_metrics.tsv
    mm_path = os.path.join(BR, 'comparisons', f'Mapping_metrics_graph_arms_{plat}', 'mapping_metrics.tsv')
    if not os.path.exists(mm_path):
        fail(f'no {mm_path}')
    prov.append((mm_path, md5(mm_path)))
    mm = {r['arm']: r for r in csv.DictReader(open(mm_path, newline=''), delimiter='\t')}
    for key in gs.MAPPING_SETUPS:
        vg_path = os.path.join(BR, f'{plat}_{gs.MAPPING_FOLDER[key]}', 'mapping', 'vg_stats_a.txt')
        if not os.path.exists(vg_path):
            fail(f'no {vg_path}')
        prov.append((vg_path, md5(vg_path)))
        text = open(vg_path).read()
        total = re.search(r'Total alignments: (\d+)', text)
        secondary = re.search(r'Total secondary: (\d+)', text)
        if not total or not secondary or int(secondary.group(1)) != 0:
            fail(f'{vg_path}: no total or secondary alignments present')
        total = int(total.group(1))
        arm = f'{plat}_{gs.CALL_FOLDER[key]}'
        if arm not in mm:
            fail(f'{mm_path}: no row {arm}')
        src_vg = f'{plat}_{gs.MAPPING_FOLDER[key]}/mapping/vg_stats_a.txt'
        for metric, pat, col in VG:
            m = re.search(pat, text)
            if not m:
                fail(f'{vg_path}: no line for {metric}')
            pct = 100 * int(m.group(1)) / total
            if abs(pct - float(mm[arm][col])) > 0.006:
                fail(f'{arm} {metric}: {pct:.4f} from vg stats, {mm[arm][col]} in mapping_metrics.tsv')
            map_rows.append([plat, key, metric, f'{pct:.4f}', src_vg])
        for metric, col in BAM:
            try:
                float(mm[arm][col])
            except ValueError:
                fail(f'{arm}: {col} is {mm[arm][col]}')
            map_rows.append([plat, key, metric, mm[arm][col], f'comparisons/Mapping_metrics_graph_arms_{plat}/mapping_metrics.tsv'])
    # BAM metrics of BWA-MEM (no graph alignment): on chr1-22, X, Y as the giraffe BAMs (bam_*), and on all contigs (bam_all_*); the same
    # reads as the giraffe set-ups, so the primary records must be equal
    arm = f'{plat}_{gs.CALL_FOLDER["bwa_grch38"]}'
    cnts = {}
    for sub, prefix in (('bam_chr1_22XY', 'bam_'), ('bam', 'bam_all_')):
        bdir = os.path.join(BR, 'comparisons', f'Mapping_metrics_graph_arms_{plat}', sub)
        bwa_path = os.path.join(bdir, f'{arm}.tsv')
        if not os.path.exists(os.path.join(bdir, f'BAM_METRICS_{arm}_DONE')) or not os.path.exists(bwa_path):
            fail(f'{arm}: no BAM_METRICS_{arm}_DONE or {bwa_path}')
        prov.append((bwa_path, md5(bwa_path)))
        cnt = cnts[prefix] = {k: int(v) for k, v in (line.rstrip('\n').split('\t') for line in open(bwa_path))}
        giraffe_primary = {int(mm[f'{plat}_{gs.CALL_FOLDER[s]}']['bam_primary']) for s in gs.MAPPING_SETUPS}
        if giraffe_primary != {cnt['primary']}:
            fail(f'{arm} ({sub}): {cnt["primary"]} primary records, the giraffe BAMs have {sorted(giraffe_primary)}')
        for metric, _ in BAM:
            if prefix == 'bam_all_' and metric == 'bam_properly_paired':
                continue
            key = metric[len('bam_'):]
            map_rows.append([plat, 'bwa_grch38', prefix + key, f'{100 * cnt[key] / cnt["primary"]:.4f}',
                             f'comparisons/Mapping_metrics_graph_arms_{plat}/{sub}/{arm}.tsv'])
    if any(cnts['bam_'][k] > cnts['bam_all_'][k] for k in ('mapped', 'properly_paired', 'mapq60')):
        fail(f'{arm}: more records on chr1-22, X, Y than on all contigs')

n_setups, n_map = len(gs.PLATFORMS) * len(gs.SETUPS), len(gs.PLATFORMS) * len(gs.MAPPING_SETUPS)
if len(rows) != n_setups * len(gs.REGIONS) * len(TYPES):
    fail(f'{len(rows)} rows, expected {n_setups * len(gs.REGIONS) * len(TYPES)}')
n_bwa = len(gs.PLATFORMS) * (2 * len(BAM) - 1)
if len(map_rows) != n_map * (len(VG) + len(BAM)) + n_bwa:
    fail(f'{len(map_rows)} alignment rows, expected {n_map * (len(VG) + len(BAM)) + n_bwa}')
outs = [write('graph_arms_data.tsv', COLS, rows), write('graph_arms_roc.tsv', ROC_COLS, roc_rows),
        write('graph_arms_mapping.tsv', MAP_COLS, map_rows)]
with open(os.path.join(HERE, 'graph_arms_data.provenance.txt'), 'w') as fh:
    for path, m in prov:
        fh.write(f'input\t{m}\t{path}\n')
    for path in (os.path.abspath(__file__), os.path.join(HERE, 'graph_setups.py')):
        fh.write(f'script\t{md5(path)}\t{path}\n')
    for out in outs:
        fh.write(f'table\t{md5(out)}\t{out}\n')
print('GRAPH_ARMS_TABLE_OK')
