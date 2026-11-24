#!/usr/bin/env python3
"""The data and the tests of Figure 3 and Supplementary Figures S8-S12: genotype and allele concordance of the nine call-set pairs
per individual (numbered 1-23) in the genome and the 10 regions, and paired Wilcoxon signed-rank tests with Holm over the 10 regions
(genome: raw p). The tests that the original run also made are checked against its posthoc_long.tsv.
Inputs: <CONCORDANCE_DIR>/summary/all_counts.tsv, <CONCORDANCE_DIR>/stats/posthoc_long.tsv. Usage: python3 concordance_table.py
"""
import csv
import hashlib
import math
import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
CD = os.environ.get('CONCORDANCE_DIR', '.')
SRC = os.path.join(CD, 'summary', 'all_counts.tsv')
POST = os.path.join(CD, 'stats', 'posthoc_long.tsv')
REGIONS = ['BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']   # CMRG and chrX PARs left out
STRATA = ['genome'] + REGIONS
TYPES = ['ALL', 'SNV', 'INDEL']
PAIRS = ['GATK_DV_illumina', 'DV_PPR_illumina', 'GATK_PPR_illumina', 'GATK_DV_aviti', 'DV_PPR_aviti', 'GATK_PPR_aviti',
         'GATK_illumina_aviti', 'DV_illumina_aviti', 'PPR_illumina_aviti']
METRIC = {'J_gt': 'genotype', 'J_al': 'allele'}
# the comparisons: x and y of each individual, as means of call-set pairs
COMPARISONS = {
    'H1_crossplatform_vs_crossmethod': ('same pipeline across platforms vs same platform across pipelines (GATK and PPR)',
                                        ['PPR_illumina_aviti', 'GATK_illumina_aviti'], ['GATK_PPR_illumina', 'GATK_PPR_aviti']),
    'H1_DV_PPR_crossplatform_vs_crossmethod': ('same pipeline across platforms vs same platform across pipelines (DV and PPR)',
                                               ['PPR_illumina_aviti', 'DV_illumina_aviti'], ['DV_PPR_illumina', 'DV_PPR_aviti']),
    'H2_PPR_vs_GATK_crossplatform': ('PPR vs GATK, concordance across platforms', ['PPR_illumina_aviti'], ['GATK_illumina_aviti']),
    'H2_DV_vs_GATK_crossplatform': ('DV vs GATK, concordance across platforms', ['DV_illumina_aviti'], ['GATK_illumina_aviti']),
    'H2_PPR_vs_DV_crossplatform': ('PPR vs DV, concordance across platforms', ['PPR_illumina_aviti'], ['DV_illumina_aviti']),
    'H3_PPRDV_vs_DVGATK_illumina': ('DV vs PPR (reference changes) vs GATK vs DV (caller changes), Illumina', ['DV_PPR_illumina'],
                                    ['GATK_DV_illumina']),
    'H3_PPRDV_vs_DVGATK_aviti': ('DV vs PPR (reference changes) vs GATK vs DV (caller changes), Aviti', ['DV_PPR_aviti'],
                                 ['GATK_DV_aviti'])}


def md5(path):
    with open(path, 'rb') as fh:
        return hashlib.md5(fh.read()).hexdigest()


def fail(msg):
    print('ERROR:', msg)
    print('CONCORDANCE_TABLE_FAILED')
    sys.exit(1)


def holm(ps):
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    adj, running = [0.0] * len(ps), 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(ps) - rank) * ps[i]))
        adj[i] = running
    return adj


def wilcoxon(x, y):
    """the test of 06_summary_stats.py: n, V (sum of the ranks of the positive differences), raw p, r = |z| / sqrt(n) at most 1,
    and the matched-pairs rank-biserial correlation"""
    d = np.asarray(x, dtype=float) - np.asarray(y, dtype=float)
    n_all = len(d)
    nz = d[d != 0]
    if len(nz) < 3:
        return n_all, float('nan'), float('nan'), float('nan'), float('nan')
    ranks = stats.rankdata(np.abs(nz))
    v = float(ranks[nz > 0].sum())
    p = float(stats.wilcoxon(nz, zero_method='wilcox', correction=True, method='approx').pvalue)
    z = stats.norm.isf(p / 2) if p > 0 else float('inf')
    rb = float((ranks[nz > 0].sum() - ranks[nz < 0].sum()) / (len(nz) * (len(nz) + 1) / 2))
    return n_all, v, p, min(abs(z) / math.sqrt(n_all), 1.0), rb


def close(a, b, tol):
    a, b = float(a), float(b)
    return (math.isnan(a) and math.isnan(b)) or abs(a - b) <= tol * max(1.0, abs(b))


# ---- concordance per individual
with open(SRC, newline='') as fh:
    src = [r for r in csv.DictReader(fh, delimiter='\t') if r['stratum'] in STRATA and r['type'] in TYPES and r['pair'] in PAIRS]
codes = sorted({r['patient'] for r in src})
if len(codes) != 23:
    fail(f'{len(codes)} individuals, expected 23')
index = {c: i + 1 for i, c in enumerate(codes)}
V = {}
data_rows = []
for r in src:
    i = index[r['patient']]
    for col, metric in METRIC.items():
        v = float(r[col])
        if not 0 <= v <= 1:
            fail(f'{col} out of range: {r["pair"]} {r["stratum"]} {r["type"]}')
        V[(metric, r['pair'], r['stratum'], r['type'], i)] = v
    data_rows.append([i, r['pair'], r['kind'], r['stratum'], r['type'], r['J_gt'], r['J_al']])
if len(data_rows) != 23 * len(PAIRS) * len(STRATA) * len(TYPES):
    fail(f'{len(data_rows)} rows, expected {23 * len(PAIRS) * len(STRATA) * len(TYPES)}')
data_rows.sort(key=lambda x: (PAIRS.index(x[1]), STRATA.index(x[3]), TYPES.index(x[4]), x[0]))

# ---- tests; those that 07_stats_figures.py of the original run also made are checked against its posthoc_long.tsv
with open(POST, newline='') as fh:
    post = {(r['comparison_id'], r['metric'], r['stratum'], r['type']): r for r in csv.DictReader(fh, delimiter='\t')
            if r['metric'] in METRIC and r['stratum'] in STRATA and r['type'] in TYPES}
unknown = sorted({k[0] for k in post} - set(COMPARISONS))
if unknown:
    fail(f'unknown comparisons in posthoc_long.tsv: {unknown}')
test_rows, checked = [], 0
for cid, (title, xs, ys) in COMPARISONS.items():
    for metric in METRIC:
        m = METRIC[metric]
        for typ in TYPES:
            res = {}
            for s in STRATA:
                x = [sum(V[(m, p, s, typ, i)] for p in xs) / len(xs) for i in range(1, 24)]
                y = [sum(V[(m, p, s, typ, i)] for p in ys) / len(ys) for i in range(1, 24)]
                n, v, p, r_eff, rb = wilcoxon(x, y)
                mx, my = 100 * float(np.mean(x)), 100 * float(np.mean(y))
                res[s] = (n, v, p, r_eff, rb, mx, my, mx - my)
                o = post.get((cid, metric, s, typ))
                if o is not None:
                    checked += 1
                    for mine, col, tol in ((mx, 'mean_x_pp', 1e-6), (my, 'mean_y_pp', 1e-6), (mx - my, 'mean_diff_pp', 1e-6),
                                           (n, 'n', 0), (v, 'V', 0), (p, 'p_raw', 1e-9), (r_eff, 'effect_r', 1e-9),
                                           (rb, 'rank_biserial', 1e-9)):
                        if not close(mine, o[col], tol):
                            fail(f'{cid} {metric} {s} {typ}: {col} {mine} here, {o[col]} in posthoc_long.tsv')
                elif cid != 'H1_DV_PPR_crossplatform_vs_crossmethod':
                    fail(f'{cid} {metric} {s} {typ}: no row in posthoc_long.tsv')
            adj = dict(zip(REGIONS, holm([res[s][2] for s in REGIONS])))
            adj['genome'] = res['genome'][2]
            for s in STRATA:
                n, v, p, r_eff, rb, mx, my, md = res[s]
                q = adj[s]
                test_rows.append([cid, title, s, typ, m, n, f'{mx:.6f}', f'{my:.6f}', f'{md:.6f}', f'{v:g}', f'{p:.6g}', f'{r_eff:.6f}',
                                  f'{rb:.6f}', 1 if s == 'genome' else 10, f'{q:.6g}', 'True' if q < 0.05 else 'False'])
print(f'{checked} tests equal to posthoc_long.tsv of 07; {len(test_rows) - checked} new')
test_rows.sort(key=lambda x: (list(COMPARISONS).index(x[0]), x[4], TYPES.index(x[3]), STRATA.index(x[2])))
expected = len(COMPARISONS) * 2 * len(TYPES) * len(STRATA)
if len(test_rows) != expected:
    fail(f'{len(test_rows)} tests, expected {expected}')

outs = []
for name, head, rows in (('concordance_data.tsv', ['individual', 'pair', 'kind', 'stratum', 'type', 'genotype_concordance',
                                                   'allele_concordance'], data_rows),
                         ('concordance_tests.tsv', ['comparison_id', 'comparison', 'stratum', 'type', 'metric', 'n', 'mean_x_pp',
                                                    'mean_y_pp', 'mean_diff_pp', 'V', 'p_raw', 'effect_r', 'rank_biserial', 'holm_m',
                                                    'p_holm', 'significant'], test_rows)):
    out = os.path.join(HERE, name)
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(head)
        w.writerows(rows)
    outs.append(out)
    print(f'{len(rows)} rows written to {out}')
with open(os.path.join(HERE, 'concordance_data.provenance.txt'), 'w') as fh:
    for path in (SRC, POST):
        fh.write(f'input\t{md5(path)}\t{path}\n')
    fh.write(f'script\t{md5(os.path.abspath(__file__))}\t{os.path.abspath(__file__)}\n')
    for out in outs:
        fh.write(f'table\t{md5(out)}\t{out}\n')
print('CONCORDANCE_TABLE_OK')
