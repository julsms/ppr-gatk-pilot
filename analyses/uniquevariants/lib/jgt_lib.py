#!/usr/bin/env python3
"""Shared parts of xplat_jgt_summary.py and xplat_jgt_novel.py: the sets, regions and types, and the test of the concordance
figures (figures/concordance/concordance_table.py): scipy wilcoxon on the non-zero differences, zero_method wilcox, continuity
correction, normal approximation; r = |z| / sqrt(n), at most 1; matched-pairs rank-biserial; Holm correction.
"""
import math

import numpy as np
from scipy import stats

SETS = ['PPRvsGATK', 'GATKvsPPR', 'PPRvsDV', 'DVvsPPR']
REGIONS = ['BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']
STRATA = ['genome'] + REGIONS
TYPES = ['ALL', 'SNV', 'INDEL']
COMPARISONS = [('PPRvsGATK', 'GATKvsPPR', 'PPR-only (vs GATK) vs GATK-only (vs PPR)'),
               ('PPRvsDV', 'DVvsPPR', 'PPR-only (vs DV) vs DV-only (vs PPR)')]
TEST_HEAD = ['comparison_id', 'comparison', 'stratum', 'type', 'n', 'mean_x_pp', 'mean_y_pp', 'mean_diff_pp', 'V', 'p_raw', 'effect_r',
             'rank_biserial', 'holm_m', 'p_holm', 'significant']


def holm(ps):
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    adj, running = [0.0] * len(ps), 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(ps) - rank) * ps[i]))
        adj[i] = running
    return adj


def wilcoxon(x, y):
    """as concordance_table.py: n, V, raw p, r = |z| / sqrt(n) at most 1, matched-pairs rank-biserial"""
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


def test_rows(fam, holm_m, p_holm):
    """one output row of the tests table: fam = [id, label, stratum, type, n, mean_x, mean_y, diff, V, p, r, rb]"""
    fmt = lambda x, f: '' if (isinstance(x, float) and math.isnan(x)) else format(x, f)
    return fam[:5] + [fmt(fam[5], '.4f'), fmt(fam[6], '.4f'), fmt(fam[7], '.4f'), fmt(fam[8], 'g'), fmt(fam[9], '.4g'),
                      fmt(fam[10], '.4f'), fmt(fam[11], '.4f'), holm_m, fmt(p_holm, '.4g'),
                      '' if math.isnan(p_holm) else p_holm < 0.05]
