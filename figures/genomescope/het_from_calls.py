#!/usr/bin/env python3
"""het_from_calls.py - heterozygous alleles per individual in the call sets, per pipeline and platform, from
uniquevariants/tables/r27_per_person.tsv of the unique-variant analysis (under PAPER_AGAIN_DIR, default .), for the Results text.
"""
import csv
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PA = os.environ.get('PAPER_AGAIN_DIR', '.')
R27 = os.path.join(PA, 'uniquevariants', 'tables', 'r27_per_person.tsv')
GENOME_KB = 2.9e6           # 2.9 Gb, the order of the GenomeScope genome length
# the side of each pipeline in a concordance pair (any pair gives the same records)
SIDE = {'GATK': ('GATK_DV', 'A'), 'DV': ('DV_PPR', 'A'), 'PPR': ('GATK_PPR', 'B')}

het = {}
with open(R27, newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['stratum'] != 'genome':
            continue
        pair, plat = r['pair'].rsplit('_', 1)
        for s, (p, side) in SIDE.items():
            if (pair, r['side']) == (p, side):
                het.setdefault((s, plat), {}).setdefault(r['patient'], 0)
                het[(s, plat)][r['patient']] += int(r['het'])
out = []
for s in ('GATK', 'DV', 'PPR'):
    for plat in ('illumina', 'aviti'):
        v = list(het[(s, plat)].values())
        if len(v) != 23:
            print('ERROR:', s, plat, len(v), 'individuals')
            print('HET_FROM_CALLS_FAILED')
            sys.exit(1)
        m = statistics.mean(v)
        out.append([s, plat, len(v), f'{m:.0f}', f'{statistics.stdev(v):.0f}', f'{m / GENOME_KB:.3f}'])
        print(f'{s} {plat}: {m / 1e6:.2f} M heterozygous alleles per individual, {m / GENOME_KB:.3f} per kb ({100 * m / GENOME_KB / 1000:.3f} %)')
with open(os.path.join(HERE, 'het_from_calls.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['call_set', 'platform', 'n', 'het_mean', 'het_sd', 'het_per_kb'])
    w.writerows(out)
print('HET_FROM_CALLS_OK')
