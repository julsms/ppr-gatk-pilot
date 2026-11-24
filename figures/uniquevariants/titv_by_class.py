#!/usr/bin/env python3
"""titv_by_class.py - Ti/Tv of the SNVs of each call set, split into the unique calls (PPR-only, GATK-only, DV-only), the shared calls,
the shared plus the unique calls of the same side (for PPR: shared + PPR-only), and all records of the side (shared,
genotype-discordant, unique and not assessed), per platform and region; also het/hom. These are the values of the Ti/Tv paragraph of
the Results (Supplementary Figure S18 draws the same ratios by region).
Method : sums of transitions and transversions over the 23 individuals (pooled ratio, as in Supplementary Figure S18), from the Ti/Tv
         and het/hom table of the unique-variant run; no new counting.
Inputs : argv[1] = <UV_DIR>/tables/r27_titv_hethom.tsv (tables.sh of analyses/uniquevariants; columns pair, platform, side, set, class,
         stratum, snv, ti, tv, het, hom)
Outputs: argv[2] = TSV (pair, platform, side, set, stratum, group, snv, titv, hethom) and the same table on stdout for the genome,
         autosomes, BB and LM; last line TITV_BY_CLASS_DONE
Usage  : python3 titv_by_class.py <UV_DIR>/tables/r27_titv_hethom.tsv titv_by_class.tsv
"""
import csv
import sys
from collections import defaultdict

GROUPS = {'unique': ['unique'], 'shared': ['shared'], 'shared+unique': ['shared', 'unique'],
          'all records': ['shared', 'gt_discordant', 'unique', 'not_assessed']}
PAIRS = ['GATK_PPR', 'DV_PPR', 'GATK_DV']
SHOW = ['genome', 'autosomes', 'BB', 'LM']

rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
acc = defaultdict(lambda: [0, 0, 0, 0, 0])
for r in rows:
    k = (r['pair'], r['platform'], r['side'], r['set'], r['stratum'], r['class'])
    acc[k] = [int(r['snv']), int(r['ti']), int(r['tv']), int(r['het']), int(r['hom'])]
keys = sorted({k[:5] for k in acc}, key=lambda k: (PAIRS.index(k[0].rsplit('_', 1)[0]), k[1] != 'illumina', k[2], k[4]))
out = []
for k in keys:
    for g, classes in GROUPS.items():
        s = [sum(acc[k + (c,)][j] for c in classes) for j in range(5)]
        titv = s[1] / s[2] if s[2] else float('nan')
        hethom = s[3] / s[4] if s[4] else float('nan')
        out.append(list(k) + [g, s[0], round(titv, 3), round(hethom, 3)])
with open(sys.argv[2], 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['pair', 'platform', 'side', 'set', 'stratum', 'group', 'snv', 'titv', 'hethom'])
    w.writerows(out)
print('pair\tplatform\tset\tstratum\tunique\tshared\tshared+unique\tall records\t(Ti/Tv; SNVs in the unique group)')
by = defaultdict(dict)
for r in out:
    by[tuple(r[:5])][r[5]] = r
for k, g in by.items():
    if k[4] in SHOW:
        print(f"{k[0]}\t{k[1]}\t{k[3]}\t{k[4]}\t{g['unique'][7]}\t{g['shared'][7]}\t{g['shared+unique'][7]}\t{g['all records'][7]}"
              f"\t{g['unique'][6]:,}")
print('TITV_BY_CLASS_DONE')
