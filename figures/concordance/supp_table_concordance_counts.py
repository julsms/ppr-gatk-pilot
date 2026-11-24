#!/usr/bin/env python3
"""The genome and BB rows of Supplementary File 4: per individual and call-set pair, the records of each call set (total, same
allele, same allele and genotype, unique, not assessed) and the two concordances, for all variants, SNVs and indels. The rows of the
other regions are those of summary/all_counts.tsv. Input: <CONCORDANCE_DIR>/summary/all_counts.tsv.
Usage: python3 supp_table_concordance_counts.py
"""
import csv
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CD = os.environ.get('CONCORDANCE_DIR', '.')
SRC = os.path.join(CD, 'summary', 'all_counts.tsv')
PAIRS = ['GATK_DV_illumina', 'DV_PPR_illumina', 'GATK_PPR_illumina', 'GATK_DV_aviti', 'DV_PPR_aviti', 'GATK_PPR_aviti',
         'GATK_illumina_aviti', 'DV_illumina_aviti', 'PPR_illumina_aviti']
STRATA = [('genome', 'genome'), ('BB', 'benchmark regions (BB)')]
TYPES = [('ALL', 'all variants'), ('SNV', 'SNVs'), ('INDEL', 'indels')]
PIPE = {'GATK': 'GATK', 'DV': 'BWA-MEM + DeepVariant', 'PPR': 'PPR'}
PLAT = {'illumina': 'Illumina', 'aviti': 'Aviti'}
COLS = ['individual', 'call_set_A', 'call_set_B', 'region', 'type', 'total_A', 'total_B', 'shared_A', 'shared_B',
        'same_genotype_A', 'same_genotype_B', 'unique_A', 'unique_B', 'not_assessed_A', 'not_assessed_B', 'allele_concordance',
        'genotype_concordance']


def md5(path):
    with open(path, 'rb') as fh:
        return hashlib.md5(fh.read()).hexdigest()


def fail(msg):
    print('ERROR:', msg)
    print('SUPP_TABLE_COUNTS_FAILED')
    sys.exit(1)


strata, types = dict(STRATA), dict(TYPES)
with open(SRC, newline='') as fh:
    src = [r for r in csv.DictReader(fh, delimiter='\t') if r['pair'] in PAIRS and r['stratum'] in strata and r['type'] in types]
codes = sorted({r['patient'] for r in src}, key=lambda c: (len(c), c))
if len(codes) != 23:
    fail(f'{len(codes)} individuals, expected 23')
rows = []
for r in src:
    for side in ('A', 'B'):
        parts = sum(int(r[f'{k}_{side}']) for k in ('GM', 'AM', 'U', 'NA'))
        if parts != int(r[f'N_{side}']):
            fail(f'{r["patient"]} {r["pair"]} {r["stratum"]} {r["type"]}: N_{side} {r[f"N_{side}"]} is not GM + AM + U + NA ({parts})')
    rows.append([r['patient'], f'{PIPE[r["a_set"]]}, {PLAT[r["a_platform"]]}', f'{PIPE[r["b_set"]]}, {PLAT[r["b_platform"]]}',
                 strata[r['stratum']], types[r['type']], r['N_A'], r['N_B'], r['S_A'], r['S_B'], r['GM_A'], r['GM_B'], r['U_A'], r['U_B'],
                 r['NA_A'], r['NA_B'], f'{float(r["J_al"]):.4f}', f'{float(r["J_gt"]):.4f}'])
expected = 23 * len(PAIRS) * len(STRATA) * len(TYPES)
if len(rows) != expected:
    fail(f'{len(rows)} rows, expected {expected}')
order = {c: i for i, c in enumerate(codes)}
rows = [x for _, x in sorted(((order[r['patient']], PAIRS.index(r['pair']), [s for s, _ in STRATA].index(r['stratum']),
                                [t for t, _ in TYPES].index(r['type'])), row) for r, row in zip(src, rows))]
out = os.path.join(HERE, 'SuppTable_concordance_counts.tsv')
with open(out, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(COLS)
    w.writerows(rows)
with open(os.path.join(HERE, 'SuppTable_concordance_counts.provenance.txt'), 'w') as fh:
    fh.write(f'input\t{md5(SRC)}\t{SRC}\n')
    fh.write(f'script\t{md5(os.path.abspath(__file__))}\t{os.path.abspath(__file__)}\n')
    fh.write(f'table\t{md5(out)}\t{out}\t{len(rows)} rows\n')
print(f'{len(rows)} rows written to {out}')
print('SUPP_TABLE_COUNTS_DONE')
