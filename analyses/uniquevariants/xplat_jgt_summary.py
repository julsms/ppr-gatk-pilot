#!/usr/bin/env python3
"""Summary and tests of the genotype concordance (J_gt) between Illumina and Aviti of the unique sets (xplat_jgt.sh). Per set,
region and type: the J_gt of the 23 individuals (mean, median, quartiles, in percent) and J_al; paired comparisons PPR vs GATK
against GATK vs PPR, and PPR vs DV against DV vs PPR, with the test of the concordance figures (lib/jgt_lib.py), Holm over the 10
regions for each comparison and type, the genome row unadjusted.
Inputs : uniquevariants/xplat_jgt/classify/<set>_illumina_aviti/<individual>.counts.tsv (UV_OUT, default uniquevariants)
Outputs: uniquevariants/xplat_jgt/tables/jgt_per_individual.tsv, jgt_summary.tsv, jgt_tests.tsv; stdout: checks and the genome
         rows; last line UNIQUE_XPLAT_JGT_SUMMARY_OK
Usage  : python3 xplat_jgt_summary.py      (in the run folder, after xplat_jgt.sh for every individual)
"""
import csv
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
from jgt_lib import COMPARISONS, REGIONS, SETS, STRATA, TEST_HEAD, TYPES, holm, test_rows, wilcoxon  # noqa: E402

N_IND = 23
O = os.path.join(os.environ.get('UV_OUT', 'uniquevariants'), 'xplat_jgt')
T = os.path.join(O, 'tables')


def fail(msg):
    print(msg)
    sys.exit('UNIQUE_XPLAT_JGT_SUMMARY_FAILED')


# per individual
rows, val, na_max = [], {}, 0.0
for s in SETS:
    files = sorted(glob.glob(os.path.join(O, 'classify', f'{s}_illumina_aviti', '*.counts.tsv')))
    if len(files) != N_IND:
        fail(f'{s}: {len(files)} counts files, expected {N_IND}')
    for f in files:
        ind = os.path.basename(f).split('.')[0]
        if not os.path.exists(f[:-len('.counts.tsv')] + '.DONE'):
            fail(f'{f}: no DONE')
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                if r['stratum'] in STRATA and r['type'] in TYPES:
                    jg, ja = float(r['J_gt']), float(r['J_al'])
                    if not (0 <= jg <= ja <= 1):
                        fail(f'{f}: J_gt {jg}, J_al {ja} out of range ({r["stratum"]} {r["type"]})')
                    na_max = max(na_max, float(r['NA_share']))
                    val[(s, r['stratum'], r['type'], ind)] = jg
                    rows.append([s, ind, r['stratum'], r['type'], r['N_A'], r['N_B'], r['GM_A'], r['AM_A'], r['U_A'], r['U_B'],
                                 r['Sbar'], r['Gbar'], r['J_al'], r['J_gt'], r['NA_share']])
inds = sorted({k[3] for k in val}, key=lambda x: int(''.join(c for c in x if c.isdigit()) or 0))
if len(inds) != N_IND:
    fail(f'{len(inds)} individuals, expected {N_IND}')
for s in SETS:
    for st in STRATA:
        for t in TYPES:
            miss = [i for i in inds if (s, st, t, i) not in val]
            if miss:
                fail(f'{s} {st} {t}: no value for {len(miss)} individuals')

# summary per set, region and type
summ = []
for s in SETS:
    for st in STRATA:
        for t in TYPES:
            x = 100 * np.array([val[(s, st, t, i)] for i in inds])
            summ.append([s, st, t, len(x), f'{x.mean():.2f}', f'{np.median(x):.2f}', f'{np.percentile(x, 25):.2f}',
                         f'{np.percentile(x, 75):.2f}'])

# tests
tests = []
for a, b, label in COMPARISONS:
    for t in TYPES:
        fam = []
        for st in STRATA:
            x = [100 * val[(a, st, t, i)] for i in inds]
            y = [100 * val[(b, st, t, i)] for i in inds]
            n, v, p, r_eff, rb = wilcoxon(x, y)
            fam.append([f'{a}_vs_{b}', label, st, t, n, np.mean(x), np.mean(y), np.mean(x) - np.mean(y), v, p, r_eff, rb])
        reg = [r for r in fam if r[2] in REGIONS]
        if len(reg) != 10:
            fail(f'{a} vs {b} {t}: {len(reg)} regions in the Holm family')
        adj = holm([r[9] for r in reg])
        for r in fam:
            if r[2] == 'genome':
                m, ph = 1, r[9]
            else:
                m, ph = 10, adj[reg.index(r)]
            tests.append(test_rows(r, m, ph))


def write(name, head, data):
    os.makedirs(T, exist_ok=True)
    with open(os.path.join(T, name), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(head)
        w.writerows(data)
    print(f'{len(data)} rows written to {name}')


write('jgt_per_individual.tsv', ['set', 'individual', 'stratum', 'type', 'N_A', 'N_B', 'GM_A', 'AM_A', 'U_A', 'U_B', 'Sbar', 'Gbar',
                                 'J_al', 'J_gt', 'NA_share'], rows)
write('jgt_summary.tsv', ['set', 'stratum', 'type', 'n', 'mean_jgt_pct', 'median_jgt_pct', 'q25_pct', 'q75_pct'], summ)
write('jgt_tests.tsv', TEST_HEAD, tests)
print(f'checks: {N_IND} individuals x {len(SETS)} sets x {len(STRATA)} strata x {len(TYPES)} types present; '
      f'J_gt <= J_al in every row; largest NA share {na_max:.5f}')
for r in summ:
    if r[1] == 'genome':
        print('summary', *r)
for r in tests:
    if r[2] == 'genome':
        print('test', *r)
print('UNIQUE_XPLAT_JGT_SUMMARY_OK')
