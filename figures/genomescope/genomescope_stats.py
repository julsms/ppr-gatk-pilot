#!/usr/bin/env python3
"""genomescope_stats.py - GenomeScope 2.0 platform comparison from the paired table of analyses/genomescope/statistical_analysis.py in
GENOMESCOPE_DIR (default genomescope): paired Wilcoxon signed-rank tests Holm-corrected over the four estimates, matched-pairs
rank-biserial r, Spearman correlations of the heterozygosity, and the table per individual (1-21), written next to this script.
"""
import csv
import hashlib
import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
GS = os.environ.get('GENOMESCOPE_DIR', 'genomescope')
SR = os.path.join(GS, 'statistical_results')
EST = [('heterozygosity', 'Heterozygosity (%)'), ('error_rate', 'Read error rate (%)'), ('model_fit', 'Model fit (%)'),
       ('genome_length_gb', 'Genome length (Gb)')]


def md5(path):
    with open(path, 'rb') as fh:
        return hashlib.md5(fh.read()).hexdigest()


def fail(msg):
    print('ERROR:', msg)
    print('GENOMESCOPE_STATS_FAILED')
    sys.exit(1)


def read(name):
    with open(os.path.join(SR, name), newline='') as fh:
        return list(csv.DictReader(fh))


def rank_biserial(x, y):
    """matched-pairs rank-biserial correlation of y - x: (sum of the ranks of the positive differences - sum of the negative ones) /
    sum of all ranks; zero differences dropped"""
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


rows = read('paired_data.csv')
old = {r['metric']: r for r in read('paired_wilcoxon_results.csv')}
data = {e: (np.array([float(r[f'illumina_{e}']) for r in rows]), np.array([float(r[f'aviti_{e}']) for r in rows])) for e, _ in EST}
n = len(rows)
if n != 21:
    fail(f'{n} pairs (expected 21 after the outlier rule)')
res = []
for e, lab in EST:
    il, av = data[e]
    w, p = stats.wilcoxon(il, av, alternative='two-sided')
    if abs(w - float(old[e]['statistic'])) > 1e-9 or abs(p - float(old[e]['p_value'])) > 1e-12 * max(1.0, p):
        fail(f'{e}: W {w} / p {p} differ from statistical_analysis.py ({old[e]["statistic"]} / {old[e]["p_value"]})')
    res.append([e, lab, n, np.median(il), *np.percentile(il, [25, 75]), np.median(av), *np.percentile(av, [25, 75]),
                int((av > il).sum()), w, p, rank_biserial(il, av), float(old[e]['effect_size_r'])])
padj = holm([r[11] for r in res])
with open(os.path.join(HERE, 'genomescope_stats.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['estimate', 'label', 'n_pairs', 'illumina_median', 'illumina_q1', 'illumina_q3', 'aviti_median', 'aviti_q1', 'aviti_q3',
                'n_aviti_higher', 'W', 'p', 'p_holm', 'rank_biserial_aviti_minus_illumina', 'statistical_analysis_r'])
    for r, pa in zip(res, padj):
        w.writerow(r[:3] + [f'{v:.4f}' for v in r[3:9]] + [r[9], f'{r[10]:g}', f'{r[11]:.3g}', f'{pa:.3g}', f'{r[12]:.3f}',
                                                            f'{r[13]:.3f}'])
        print(f'{r[1]}: Illumina {r[3]:.3f} [{r[4]:.3f}-{r[5]:.3f}], Aviti {r[6]:.3f} [{r[7]:.3f}-{r[8]:.3f}], Aviti higher in {r[9]}/{n}, '
              f'W {r[10]:g}, p {r[11]:.2g}, Holm {pa:.2g}, rank-biserial {r[12]:.3f} (r of statistical_analysis.py {r[13]:.3f})')
cor = []
for pi, plat in ((0, 'Illumina'), (1, 'Aviti')):
    for e, _ in EST[1:]:
        rho, p = stats.spearmanr(data['heterozygosity'][pi], data[e][pi])
        cor.append([plat, 'heterozygosity', e, n, f'{rho:.3f}', f'{p:.3g}'])
        print(f'{plat}: heterozygosity vs {e}: Spearman rho {rho:.3f}, p {p:.2g}')
d_het = data['heterozygosity'][1] - data['heterozygosity'][0]
d_err = data['error_rate'][1] - data['error_rate'][0]
rho, p = stats.spearmanr(d_het, d_err)
cor.append(['Aviti minus Illumina', 'heterozygosity difference', 'error_rate difference', n, f'{rho:.3f}', f'{p:.3g}'])
print(f'Aviti minus Illumina: heterozygosity difference vs error rate difference: Spearman rho {rho:.3f}, p {p:.2g}')
with open(os.path.join(HERE, 'genomescope_correlations.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['platform', 'x', 'y', 'n', 'spearman_rho', 'p'])
    w.writerows(cor)
with open(os.path.join(HERE, 'genomescope_per_individual.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['individual', 'platform'] + [e for e, _ in EST])
    for i in range(n):
        for pi, plat in ((0, 'illumina'), (1, 'aviti')):
            w.writerow([i + 1, plat] + [f'{data[e][pi][i]:.6g}' for e, _ in EST])
outl = read('outliers_removed.csv')
print(f'outliers removed by the IQR x 3 rule of statistical_analysis.py (all four estimates, both platforms pooled): {len(outl)} rows, '
      f'{len({r["sample"].replace("ICDG", "") for r in outl})} individuals')
with open(os.path.join(HERE, 'genomescope_stats.provenance.txt'), 'w') as fh:
    for name in ('paired_data.csv', 'paired_wilcoxon_results.csv', 'outliers_removed.csv'):
        fh.write(f'input\t{md5(os.path.join(SR, name))}\t{os.path.join(SR, name)}\n')
    fh.write(f'script\t{md5(os.path.abspath(__file__))}\t{os.path.abspath(__file__)}\n')
    for name in ('genomescope_stats.tsv', 'genomescope_correlations.tsv', 'genomescope_per_individual.tsv'):
        fh.write(f'table\t{md5(os.path.join(HERE, name))}\t{os.path.join(HERE, name)}\n')
print('GENOMESCOPE_STATS_OK')
