#!/usr/bin/env python3
"""Genotype concordance (J_gt) between Illumina and Aviti of the novel variants of each unique set, from the classes of xplat_jgt.sh:
the formula of 04_classify_count.py, J_gt = Gbar / (Sbar + U_A + U_B), with S = GM + AM and G = GM of each side and Sbar, Gbar their
means over the two sides, on the records of each side that are novel (autosomes; absent from dbSNP b157, gnomAD genomes and
1000 Genomes: the novel column of vep/vep_sites.tsv.gz). A record is in a region when its key table flag is 1; types ALL, SNV, INDEL.
Self-check: the same counting on all records must give the J_gt of 04 (counts.tsv) in every row; the script stops otherwise.
Tests: as xplat_jgt_summary.py (lib/jgt_lib.py), on the individuals with a defined J_gt in both sets (a region with no novel record
has none).
Inputs : uniquevariants/xplat_jgt/classify/<set>_illumina_aviti/<individual>.{records.tsv.gz,counts.tsv}, uniquevariants/xplat_jgt/
         prep/<set>_<platform>/<individual>.keys.tsv.gz, uniquevariants/vep/vep_sites.tsv.gz (UV_OUT, default uniquevariants)
Outputs: uniquevariants/xplat_jgt/tables/novel_jgt_per_individual.tsv, novel_jgt_summary.tsv, novel_jgt_tests.tsv; stdout: checks
         and the genome rows; last line NOVEL_JGT_DONE
Usage  : python3 xplat_jgt_novel.py        (in the run folder, after xplat_jgt.sh for every individual; one core, a few minutes)
"""
import csv
import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
from jgt_lib import COMPARISONS, REGIONS, SETS, STRATA, TEST_HEAD, TYPES, holm, test_rows, wilcoxon  # noqa: E402

UV = os.environ.get('UV_OUT', 'uniquevariants')
O = os.path.join(UV, 'xplat_jgt')
T = os.path.join(O, 'tables')
VEP = os.path.join(UV, 'vep', 'vep_sites.tsv.gz')
TOL = 1e-6


def fail(msg):
    print(msg)
    sys.exit('NOVEL_JGT_FAILED')


def keyed(df):
    df['key'] = df['CHROM'] + ':' + df['POS'] + ':' + df['REF'] + ':' + df['ALT']
    return df


def side_table(rec, keys, side):
    """records of one side with class, type, AUTO and the region flags of its key table"""
    s = rec[rec['side'] == side][['key', 'class']]
    k = keys.set_index('key')
    missing = ~s['key'].isin(k.index)
    if missing.any():
        fail(f'side {side}: {int(missing.sum())} record keys not in the key table')
    return s.join(k[['TYPE', 'AUTO'] + REGIONS], on='key')


def counts(sa, sb):
    """per stratum and type: Sbar, Gbar, U_A, U_B and J_gt as 04 (NA records left out of the denominator)"""
    out = {}
    for st in STRATA:
        for t in TYPES:
            c = {}
            for nm, s in (('A', sa), ('B', sb)):
                m = s if st == 'genome' else s[s[st] == '1']
                if t != 'ALL':
                    m = m[m['TYPE'] == t]
                v = m['class'].value_counts()
                c[nm] = {x: int(v.get(x, 0)) for x in ('GM', 'AM', 'U', 'NA')}
            sbar = (c['A']['GM'] + c['A']['AM'] + c['B']['GM'] + c['B']['AM']) / 2
            gbar = (c['A']['GM'] + c['B']['GM']) / 2
            den = sbar + c['A']['U'] + c['B']['U']
            out[(st, t)] = (sbar, gbar, c['A']['U'], c['B']['U'], gbar / den if den > 0 else float('nan'), c['A'], c['B'])
    return out


vep = pd.read_csv(VEP, sep='\t', dtype=str, keep_default_na=False, usecols=['key', 'novel'])
novel = set(vep.loc[vep['novel'] == '1', 'key'])
print(f'novel alleles in vep_sites: {len(novel):,}')
rows, val, n_check, max_diff = [], {}, 0, 0.0
for s in SETS:
    pair = f'{s}_illumina_aviti'
    cdir = os.path.join(O, 'classify', pair)
    inds = sorted(f.split('.')[0] for f in os.listdir(cdir) if f.endswith('.records.tsv.gz'))
    for ind in inds:
        rec = keyed(pd.read_csv(os.path.join(cdir, f'{ind}.records.tsv.gz'), sep='\t', dtype=str, keep_default_na=False,
                                usecols=['side', 'CHROM', 'POS', 'REF', 'ALT', 'class']))
        kt = {p: keyed(pd.read_csv(os.path.join(O, 'prep', f'{s}_{p}', f'{ind}.keys.tsv.gz'), sep='\t', dtype=str,
                                   keep_default_na=False, usecols=['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'AUTO'] + REGIONS))
              for p in ('illumina', 'aviti')}
        sa, sb = side_table(rec, kt['illumina'], 'A'), side_table(rec, kt['aviti'], 'B')
        # self-check against 04 on all records
        ref = pd.read_csv(os.path.join(cdir, f'{ind}.counts.tsv'), sep='\t', dtype=str, keep_default_na=False)
        # 04 writes NA where J_gt is undefined (no record in the stratum and type); it must be undefined here too
        ref = {(r['stratum'], r['type']): (float('nan') if r['J_gt'] == 'NA' else float(r['J_gt'])) for _, r in ref.iterrows()}
        for k, v in counts(sa, sb).items():
            if math.isnan(v[4]) or math.isnan(ref[k]):
                if not (math.isnan(v[4]) and math.isnan(ref[k])):
                    fail(f'{pair} {ind} {k}: J_gt {v[4]} here, {ref[k]} in 04: only one is undefined')
            else:
                d = abs(v[4] - ref[k])
                max_diff = max(max_diff, d)
                if d > TOL:
                    fail(f'{pair} {ind} {k}: J_gt {v[4]:.6f} here, {ref[k]:.6f} in 04: the counting differs')
            n_check += 1
        na = sa[sa['key'].isin(novel) & (sa['AUTO'] == '1')]
        nb = sb[sb['key'].isin(novel) & (sb['AUTO'] == '1')]
        for (st, t), v in counts(na, nb).items():
            val[(s, st, t, ind)] = v[4]
            rows.append([s, ind, st, t, v[5]['GM'] + v[5]['AM'] + v[5]['U'] + v[5]['NA'],
                         v[6]['GM'] + v[6]['AM'] + v[6]['U'] + v[6]['NA'], v[5]['GM'], v[5]['AM'], v[2], v[3], f'{v[0]:g}',
                         f'{v[1]:g}', '' if math.isnan(v[4]) else f'{v[4]:.6f}'])
print(f'self-check: J_gt of 04 reproduced on all records in {n_check} rows (largest difference {max_diff:.2e})')

inds = sorted({k[3] for k in val})
summ = []
for s in SETS:
    for st in STRATA:
        for t in TYPES:
            x = 100 * np.array([val[(s, st, t, i)] for i in inds if not math.isnan(val[(s, st, t, i)])])
            summ.append([s, st, t, len(x)] + ([f'{x.mean():.2f}', f'{np.median(x):.2f}', f'{np.percentile(x, 25):.2f}',
                                               f'{np.percentile(x, 75):.2f}'] if len(x) else ['', '', '', '']))
tests = []
for a, b, label in COMPARISONS:
    for t in TYPES:
        fam = []
        for st in STRATA:
            both = [i for i in inds if not (math.isnan(val[(a, st, t, i)]) or math.isnan(val[(b, st, t, i)]))]
            x = [100 * val[(a, st, t, i)] for i in both]
            y = [100 * val[(b, st, t, i)] for i in both]
            n, v, p, r_eff, rb = wilcoxon(x, y) if len(both) >= 3 else (len(both),) + (float('nan'),) * 4
            fam.append([f'{a}_vs_{b}', label, st, t, n, np.mean(x) if x else float('nan'), np.mean(y) if y else float('nan'),
                        (np.mean(x) - np.mean(y)) if x else float('nan'), v, p, r_eff, rb])
        reg = [r for r in fam if r[2] in REGIONS and not math.isnan(r[9])]
        adj = dict(zip([r[2] for r in reg], holm([r[9] for r in reg]))) if reg else {}
        for r in fam:
            m, ph = (1, r[9]) if r[2] == 'genome' else (len(reg), adj.get(r[2], float('nan')))
            tests.append(test_rows(r, m, ph))


def write(name, head, data):
    os.makedirs(T, exist_ok=True)
    with open(os.path.join(T, name), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(head)
        w.writerows(data)
    print(f'{len(data)} rows written to {name}')


write('novel_jgt_per_individual.tsv', ['set', 'individual', 'stratum', 'type', 'N_A', 'N_B', 'GM_A', 'AM_A', 'U_A', 'U_B', 'Sbar',
                                       'Gbar', 'J_gt'], rows)
write('novel_jgt_summary.tsv', ['set', 'stratum', 'type', 'n_defined', 'mean_jgt_pct', 'median_jgt_pct', 'q25_pct', 'q75_pct'], summ)
write('novel_jgt_tests.tsv', TEST_HEAD, tests)
for r in summ:
    if r[1] == 'genome':
        print('summary', *r)
for r in tests:
    if r[2] == 'genome':
        print('test', *r)
print('NOVEL_JGT_DONE')
