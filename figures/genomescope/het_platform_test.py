#!/usr/bin/env python3
"""het_platform_test.py - does the heterozygosity of the variant calls differ between Illumina and Aviti in the same individuals? Set
against the GenomeScope heterozygosity, which was higher on Aviti in every individual. Per individual and pipeline: the heterozygous
rate per kb of a 2.9 Gb genome (as in het_from_calls.py) and the het/hom ratio, genome stratum, SNVs and small indels. Paired Wilcoxon
signed-rank test Illumina vs Aviti, Holm-corrected over the three pipelines for each measure, matched-pairs rank-biserial correlation
as the effect size (as in genomescope_stats.py); on all 23 individuals and on the 21 of the GenomeScope comparison (individuals 1 and
24 left out there: their fits failed on both platforms). The Results give the values of the 21.
Inputs : <PAPER_AGAIN_DIR>/uniquevariants/tables/r27_per_person.tsv (the counts of het_from_calls.py; default PAPER_AGAIN_DIR = .)
Outputs: het_platform_test.tsv (set of individuals, pipeline, measure, n, Illumina and Aviti medians, individuals higher on Aviti,
         median Aviti minus Illumina, W, p, Holm p, rank-biserial r) in this folder; last line HET_PLATFORM_TEST_OK
Usage  : python3 het_platform_test.py
"""
import csv
import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
PA = os.environ.get('PAPER_AGAIN_DIR', '.')
R27 = os.path.join(PA, 'uniquevariants', 'tables', 'r27_per_person.tsv')
GENOME_KB = 2.9e6
SIDE = {'GATK': ('GATK_DV', 'A'), 'DV': ('DV_PPR', 'A'), 'PPR': ('GATK_PPR', 'B')}   # as in het_from_calls.py
GS_OUT = {'1', '24'}                                                             # GenomeScope failed fits (individual numbers)


def rank_biserial(x, y):
    d = np.asarray(y, float) - np.asarray(x, float)
    d = d[d != 0]
    r = stats.rankdata(np.abs(d))
    return (r[d > 0].sum() - r[d < 0].sum()) / r.sum()


def holm(ps):
    order, adj, run = sorted(range(len(ps)), key=lambda i: ps[i]), [0.0] * len(ps), 0.0
    for k, i in enumerate(order):
        run = max(run, (len(ps) - k) * ps[i])
        adj[i] = min(1.0, run)
    return adj


cnt = {}
with open(R27, newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['stratum'] != 'genome':
            continue
        pair, plat = r['pair'].rsplit('_', 1)
        for s, (p, side) in SIDE.items():
            if (pair, r['side']) == (p, side):
                c = cnt.setdefault((s, plat, r['patient']), [0, 0])
                c[0] += int(r['het'])
                c[1] += int(r['hom'])
people = sorted({k[2] for k in cnt})
if len(people) != 23:
    print(f'ERROR: {len(people)} individuals, expected 23')
    print('HET_PLATFORM_TEST_FAILED')
    sys.exit(1)
rows = []
for subset, keep in (('all 23', people), ('GenomeScope 21', [p for p in people if p.replace('ICDG', '') not in GS_OUT])):
    for measure in ('het_per_kb', 'het_hom'):
        block = []
        for s in ('GATK', 'DV', 'PPR'):
            val = {}
            for plat in ('illumina', 'aviti'):
                h = np.array([cnt[(s, plat, p)][0] for p in keep], float)
                o = np.array([cnt[(s, plat, p)][1] for p in keep], float)
                val[plat] = h / GENOME_KB if measure == 'het_per_kb' else h / o
            il, av = val['illumina'], val['aviti']
            w, p = stats.wilcoxon(il, av, alternative='two-sided')
            block.append([subset, s, measure, len(keep), np.median(il), np.median(av), int((av > il).sum()), np.median(av - il), w, p,
                          rank_biserial(il, av)])
        for b, pa in zip(block, holm([b[9] for b in block])):
            rows.append(b + [pa])
with open(os.path.join(HERE, 'het_platform_test.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['individuals', 'pipeline', 'measure', 'n', 'illumina_median', 'aviti_median', 'n_aviti_higher', 'median_aviti_minus_illumina',
                'W', 'p', 'p_holm', 'rank_biserial_aviti_minus_illumina'])
    for r in rows:
        w.writerow(r[:4] + [f'{v:.4f}' for v in r[4:6]] + [r[6], f'{r[7]:.4f}', f'{r[8]:g}', f'{r[9]:.3g}', f'{r[11]:.3g}', f'{r[10]:.3f}'])
        print(f'{r[0]} {r[1]} {r[2]}: Illumina {r[4]:.4f}, Aviti {r[5]:.4f}, Aviti higher in {r[6]}/{r[3]}, median difference {r[7]:.4f}, '
              f'W {r[8]:g}, p {r[9]:.3g}, Holm {r[11]:.3g}, rank-biserial {r[10]:.3f}')
print('HET_PLATFORM_TEST_OK')
