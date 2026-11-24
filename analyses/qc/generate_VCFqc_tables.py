#!/usr/bin/env python3
"""VCF QC per individual of the six call sets (GATK, DV and PPR on Illumina and Aviti), from the variant alleles of the concordance
analysis (PASS, multi-allelic records split, normalized, carried by the individual): alleles, SNVs and indels of the genome stratum
(concordance/summary/all_counts.tsv), Ti/Tv and het/hom (uniquevariants/tables/r27_per_person.tsv), and singletons, the alleles that
no other of the 23 individuals carries in the same call set (the key tables in concordance/prep/<set>_<platform>/).
Run in the QC folder; PAPER_AGAIN_DIR (default ..) holds concordance/ and uniquevariants/. Output: tables_vcf_qc/.
"""
import csv
import os
import re
import statistics

import numpy as np
import pandas as pd

PA = os.environ.get('PAPER_AGAIN_DIR', '..')
COUNTS = os.path.join(PA, 'concordance', 'summary', 'all_counts.tsv')
R27 = os.path.join(PA, 'uniquevariants', 'tables', 'r27_per_person.tsv')
PREP = os.path.join(PA, 'concordance', 'prep')
OUTPUT_DIR = 'tables_vcf_qc'
SETS = ['GATK', 'DV', 'PPR']
PLATFORMS = ['illumina', 'aviti']
RELATED = '6'      # one of the two first-degree relatives: left out of the singleton summary
# the side of each call set in the pairs of the concordance run, as in figures/qc/table1_revised.py
SIDES = {'GATK': [('GATK_DV', 'A'), ('GATK_PPR', 'A')], 'DV': [('GATK_DV', 'B'), ('DV_PPR', 'A')],
         'PPR': [('DV_PPR', 'B'), ('GATK_PPR', 'B')]}
TYPES = {'ALL': 'total_variants', 'SNV': 'snv_count', 'INDEL': 'indel_count'}
COLUMNS = ['total_variants', 'snv_count', 'indel_count', 'titv_ratio', 'het_hom_ratio', 'singleton_count']


def individual(patient):
    """the number of the individual in a patient name of the concordance tables"""
    return re.search(r'\d+', patient).group()


def totals():
    """(set, platform, patient) -> alleles, SNVs and indels of the genome stratum; a call set has the same totals in every pair"""
    tot = {}
    with open(COUNTS, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['stratum'] != 'genome' or r['type'] not in TYPES:
                continue
            pair, plat = r['pair'].rsplit('_', 1)
            for s, sides in SIDES.items():
                for p, side in sides:
                    if p == pair and plat in PLATFORMS:
                        d = tot.setdefault((s, plat, r['patient']), {})
                        v = int(r[f'N_{side}'])
                        if d.setdefault(TYPES[r['type']], v) != v:
                            raise SystemExit(f'{s} {plat}: the totals of one individual differ between pairs')
    return tot


def classes():
    """(set, platform, patient) -> alleles, Ti, Tv, het and hom: the four concordance classes of one side summed"""
    acc = {}
    with open(R27, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['stratum'] != 'genome':
                continue
            pair, plat = r['pair'].rsplit('_', 1)
            for s, sides in SIDES.items():
                if (pair, r['side']) == sides[0]:
                    a = acc.setdefault((s, plat, r['patient']), {'n': 0, 'ti': 0, 'tv': 0, 'het': 0, 'hom': 0})
                    for c in a:
                        a[c] += int(r[c])
    return acc


def singletons(s, plat, patients):
    """per patient: the alleles (CHROM, POS, REF, ALT) of its key table that no other patient of the list has, the rows of its key
    table, and the alleles that occur twice in it"""
    keys, rows = {}, {}
    for p in patients:
        k = pd.read_csv(os.path.join(PREP, f'{s}_{plat}', f'{p}.keys.tsv.gz'), sep='\t', dtype=str, keep_default_na=False,
                        usecols=['CHROM', 'POS', 'REF', 'ALT'])
        rows[p] = len(k)
        keys[p] = np.unique(pd.util.hash_pandas_object(k, index=False).to_numpy())     # one 64-bit hash per allele
    alleles, carriers = np.unique(np.concatenate(list(keys.values())), return_counts=True)
    once = alleles[carriers == 1]
    single = {p: int(np.isin(h, once, assume_unique=True).sum()) for p, h in keys.items()}
    return single, rows, {p: rows[p] - len(h) for p, h in keys.items()}


def summary_row(rows, col, skip=()):
    v = [r[col] for r in rows if r['sample'] not in skip]
    return [col, round(statistics.mean(v), 3), round(statistics.stdev(v), 3), round(min(v), 3), round(max(v), 3),
            round(statistics.median(v), 3), len(v)]


tot, cls = totals(), classes()
os.makedirs(OUTPUT_DIR, exist_ok=True)
for s in SETS:
    for plat in PLATFORMS:
        patients = sorted({k[2] for k in tot if k[:2] == (s, plat)}, key=lambda p: int(individual(p)))
        if len(patients) != 23:
            raise SystemExit(f'{s} {plat}: {len(patients)} individuals in all_counts.tsv')
        single, n_rows, n_twice = singletons(s, plat, patients)
        rows = []
        for p in patients:
            t, a = tot[(s, plat, p)], cls[(s, plat, p)]
            if not a['n'] == n_rows[p] == t['total_variants']:
                raise SystemExit(f'{s} {plat}: the alleles of one individual differ between the input tables')
            rows.append({'sample': individual(p), **t, 'titv_ratio': a['ti'] / a['tv'], 'het_hom_ratio': a['het'] / a['hom'],
                         'singleton_count': single[p]})
        name = f'{plat}_{s.lower()}'
        with open(os.path.join(OUTPUT_DIR, f'vcf_qc_per_sample_{name}.tsv'), 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t', lineterminator='\n')
            w.writerow(['sample'] + COLUMNS)
            for r in rows:
                w.writerow([r['sample']] + [round(r[c], 4) if c.endswith('ratio') else r[c] for c in COLUMNS])
        with open(os.path.join(OUTPUT_DIR, f'vcf_qc_summary_{name}.tsv'), 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t', lineterminator='\n')
            w.writerow(['metric', 'mean', 'sd', 'min', 'max', 'median', 'n', 'note'])
            for c in COLUMNS:
                if c == 'singleton_count':
                    w.writerow(summary_row(rows, c, skip={RELATED}) + [f'individual {RELATED} (a first-degree relative) left out'])
                else:
                    w.writerow(summary_row(rows, c) + [''])
        print(f'{s} {plat}: {len(rows)} individuals, {sum(n_twice.values())} alleles twice in a key table -> {OUTPUT_DIR}/*_{name}.tsv')
