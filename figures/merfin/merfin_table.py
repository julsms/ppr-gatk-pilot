#!/usr/bin/env python3
"""merfin_table.py - data of Figure 5 and of the Merfin supplementary figures, from the outputs of analyses/kmers in KV_DIR (default
kmers_validation/revision): mf_rates.tsv (cohort, direct mode, individuals numbered 1-23), mf_hg002.tsv, mf_hg002_regions.tsv and
mf_data.provenance.txt, written next to this script. Stops when a check fails (all 48 samples done, 23 individuals in every cell).
"""
import csv
import hashlib
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
KV = os.environ.get('KV_DIR', 'kmers_validation/revision')
PLATFORMS = ['illumina', 'aviti']
REGIONS = ['BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']   # the order of Figure 1
STRATA = ['genome'] + REGIONS
TYPES = ['ALL', 'SNV', 'INDEL']
SETS = ['GATK_ALL', 'DV_ALL', 'PPR_ALL', 'GATKvsPPR', 'PPRvsGATK', 'DVvsPPR', 'PPRvsDV', 'PPRvsBOTH']
HG002_FULL = {'GATK_ALL': 'all_GATK', 'DV_ALL': 'all_DVH002', 'PPR_ALL': 'all_PPRH002'}   # the other sets: <set>_<platform>
CLASSES = ['TP', 'FP', 'FP_gt', 'NE']


def md5(path):
    with open(path, 'rb') as fh:
        return hashlib.md5(fh.read()).hexdigest()


def fail(msg):
    print('ERROR:', msg)
    print('MF_TABLE_FAILED')
    sys.exit(1)


inputs = []


def need(path):
    if not os.path.isfile(path):
        fail(f'missing input {path}')
    inputs.append(path)
    return path


def read(path):
    with open(need(path), newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write(name, head, rows):
    path = os.path.join(HERE, name)
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(head)
        w.writerows(rows)
    print(f'{len(rows)} rows written to {path}')
    return path


def f4(v):
    return f'{v:.4f}'


for p in [os.path.join(KV, 'tables', 'TABLES_DONE'), os.path.join(KV, 'hg002', 'HG002_MERFIN_DONE')] + [
        os.path.join(KV, 'hg002', 'kstar', p, 'KSTAR_DONE') for p in PLATFORMS]:
    need(p)

# ---- every sample done
comp = read(os.path.join(KV, 'tables', 'completeness.tsv'))
open_rows = [f"{r['platform']} {r['kind']}" for r in comp if r['sample_done'] != 'yes']
if open_rows or len(comp) != 2 * 24:
    fail(f'completeness: {len(comp)} rows, not done: {open_rows}')

# ---- cohort: the Merfin counts per individual (direct mode)
C = {}
for r in read(os.path.join(KV, 'tables', 'counts_long.tsv')):
    if (r['sample_kind'] != 'cohort' or r['mode'] != 'direct' or r['set_id'] not in SETS or r['stratum'] not in STRATA
            or r['type'] not in TYPES):
        continue
    n, v = int(r['n']), int(r['validated'])
    if v > n:
        fail(f'{r}: validated above n')
    if n and abs(v / n - float(r['rate'])) > 1e-5:
        fail(f'{r}: the rate column differs from validated / n')
    key = (r['platform'], r['set_id'], r['stratum'], r['type'], r['sample'])
    if key in C:
        fail(f'{key}: two rows')
    C[key] = (n, v)
codes = sorted({k[4] for k in C}, key=lambda c: (len(c), c))
if len(codes) != 23:
    fail(f'{len(codes)} individuals, expected 23')
index = {c: i + 1 for i, c in enumerate(codes)}
rates = []
for plat in PLATFORMS:
    for s in SETS:
        for st in STRATA:
            for t in TYPES:
                for c in codes:
                    if (plat, s, st, t, c) not in C:
                        fail(f'{plat} {s} {st} {t}: individual {index[c]} has no counts')
                    n, v = C[(plat, s, st, t, c)]
                    rates.append([plat, s, st, t, index[c], n, v])

# ---- HG002: Merfin x GIAB, the k-mer context, and the regions
M2 = defaultdict(int)
for r in read(os.path.join(KV, 'hg002', 'merfin_giab_2x2.tsv')):
    M2[(r['uset'], r['platform'], r['mode'], r['merfin'], r['giab'])] += int(r['n'])
SUMM = {(r['uset'], r['platform'], r['mode']): r for r in read(os.path.join(KV, 'hg002', 'merfin_giab_summary.tsv'))}
K2, KS = defaultdict(int), {}
for plat in PLATFORMS:
    for r in read(os.path.join(KV, 'hg002', 'kstar', plat, 'kstar_2x2.tsv')):
        if r['stratum'] == 'genome' and r['type'] == 'ALL':
            K2[(r['uset'], r['platform'], r['mode'], r['flag_name'], r['flag'], r['merfin'], r['giab'])] += int(r['n'])
    for r in read(os.path.join(KV, 'hg002', 'kstar', plat, 'kstar_summary.tsv')):
        KS[(r['uset'], r['platform'], r['mode'], r['stratum'], r['type'])] = r
hg, regions = [], []
for plat in PLATFORMS:
    for s in SETS:
        u = HG002_FULL.get(s, f'{s}_{plat}')
        md = 'context' if s in HG002_FULL else 'direct'
        cell = {(m, g): M2[(u, plat, md, m, g)] for m in '01' for g in CLASSES}
        tot = sum(cell.values())
        ev = tot - cell[('0', 'NE')] - cell[('1', 'NE')]
        if not ev:
            fail(f'HG002 {s} {plat}: no record inside the benchmark')
        prec = 100 * (cell[('0', 'TP')] + cell[('1', 'TP')]) / ev
        merf = 100 * sum(cell[('1', g)] for g in ('TP', 'FP', 'FP_gt')) / ev
        for src, row in (('merfin_giab_summary', SUMM.get((u, plat, md))), ('kstar_summary', KS.get((u, plat, md, 'genome', 'ALL')))):
            if row is None:
                fail(f'HG002 {s} {plat}: no row in {src}')
            if (int(row['records']) != tot or int(row['evaluable']) != ev or abs(prec - 100 * float(row['giab_precision'])) > 0.01
                    or abs(merf - 100 * float(row['merfin_share_evaluable'])) > 0.01):
                fail(f'HG002 {s} {plat}: the 2x2 ({tot}, {ev}, {prec:.3f}, {merf:.3f}) differs from {src}')
        sup = []
        for g in ('TP', 'FP', 'FP_gt'):
            k = cell[('0', g)] + cell[('1', g)]
            sup.append(f4(100 * cell[('1', g)] / k) if k else 'NA')
        ctx = []
        for gl in (('TP',), ('FP', 'FP_gt')):
            total = sum(cell[('1', x)] for x in gl)
            ks = sum(K2[(u, plat, md, 'kstar', '1', '1', x)] for x in gl)
            rep = sum(K2[(u, plat, md, 'paralog_like', '1', '1', x)] for x in gl)
            covered = sum(K2[(u, plat, md, 'paralog_like', f, '1', x)] for f in '01' for x in gl)
            if covered != total or ks > rep:
                fail(f'HG002 {s} {plat} {gl}: the k-mer flags do not add up ({covered} of {total}; kstar {ks}, repeated {rep})')
            ctx += [f4(100 * ks / total), f4(100 * (rep - ks) / total), f4(100 * (total - rep) / total)] if total else ['NA'] * 3
        hg.append([plat, s, tot, ev, f4(prec), f4(merf), f4(merf - prec),
                   f4(100 * (cell[('0', 'NE')] + cell[('1', 'NE')]) / tot)] + sup + ctx)
        for st in STRATA:
            row = KS.get((u, plat, md, st, 'ALL'))
            if row is None:
                fail(f'HG002 {s} {plat} {st}: no row in kstar_summary')
            rec, ev2 = int(row['records']), int(row['evaluable'])
            regions.append([plat, s, st, rec, ev2, f4(100 * float(row['giab_precision'])) if ev2 else 'NA',
                            f4(100 * float(row['merfin_share_evaluable'])) if ev2 else 'NA', f4(100 * ev2 / rec) if rec else 'NA'])

outs = [write('mf_rates.tsv', ['platform', 'set', 'region', 'type', 'individual', 'n', 'validated'], rates),
        write('mf_hg002.tsv', ['platform', 'set', 'records', 'evaluable', 'precision', 'merfin', 'gap', 'outside', 'sup_TP', 'sup_FP',
                               'sup_FP_gt', 'ctx_TP_kstar', 'ctx_TP_repeat', 'ctx_TP_single', 'ctx_FP_kstar', 'ctx_FP_repeat',
                               'ctx_FP_single'], hg),
        write('mf_hg002_regions.tsv', ['platform', 'set', 'region', 'records', 'evaluable', 'precision', 'merfin', 'inside'], regions)]
with open(os.path.join(HERE, 'mf_data.provenance.txt'), 'w') as fh:
    for path in inputs:
        fh.write(f'input\t{md5(path)}\t{path}\n')
    fh.write(f'script\t{md5(os.path.abspath(__file__))}\t{os.path.abspath(__file__)}\n')
    for path in outs:
        fh.write(f'table\t{md5(path)}\t{path}\n')
print('MF_TABLE_OK')
