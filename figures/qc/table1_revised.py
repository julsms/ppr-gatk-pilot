#!/usr/bin/env python3
"""Supplementary Table S2: the FASTQ and BAM rows of the QC table (analyses/qc, Table1_manuscript.tsv; DV used the GATK alignments)
and the VCF rows of GATK, DV and PPR from the variant alleles of the concordance analysis (PASS, split, normalized, carried by the
individual), as mean ± SD (range) over the 23 individuals, variant counts in millions.
Inputs: <PAPER_AGAIN_DIR>/qc/Table1_manuscript.tsv, concordance/summary/all_counts.tsv, uniquevariants/tables/r27_per_person.tsv.
Output: Table1_revised.tsv and Table1_revised.provenance.txt (md5 of the inputs) next to this script.
"""
import csv
import hashlib
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PA = os.environ.get('PAPER_AGAIN_DIR', '.')
QC = os.path.join(PA, 'qc', 'Table1_manuscript.tsv')
COUNTS = os.path.join(PA, 'concordance', 'summary', 'all_counts.tsv')
R27 = os.path.join(PA, 'uniquevariants', 'tables', 'r27_per_person.tsv')
SETS = ['GATK', 'DV', 'PPR']
PLATFORMS = ['illumina', 'aviti']
COLS = [f'{pl.capitalize()} {s}' for pl in PLATFORMS for s in SETS]
# the side of each call set in a pair of the concordance run, to read its totals and its Ti/Tv, het/hom
SIDES = {'GATK': [('GATK_DV', 'A'), ('GATK_PPR', 'A')], 'DV': [('GATK_DV', 'B'), ('DV_PPR', 'A')],
         'PPR': [('DV_PPR', 'B'), ('GATK_PPR', 'B')]}


def md5(path):
    with open(path, 'rb') as fh:
        return hashlib.md5(fh.read()).hexdigest()


def fail(msg):
    print('ERROR:', msg)
    print('TABLE1_REVISED_FAILED')
    sys.exit(1)


def summary(v, dec=0, ratio=False, scale=1.0):
    v = [x / scale for x in v]
    m, sd = statistics.mean(v), statistics.stdev(v)
    if ratio:
        return f'{m:.2f} ± {sd:.2f}'
    return f'{m:,.{dec}f} ± {sd:,.{dec}f} ({min(v):,.{dec}f}-{max(v):,.{dec}f})'


# ---- totals per person (genome; all variants, SNVs, indels)
tot = {}
with open(COUNTS, newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['stratum'] != 'genome' or r['type'] not in ('ALL', 'SNV', 'INDEL'):
            continue
        pair, plat = r['pair'].rsplit('_', 1)
        for s, sides in SIDES.items():
            for p, side in sides:
                if p == pair and plat in PLATFORMS:
                    k = (s, plat, r['type'], r['patient'])
                    v = int(r[f'N_{side}'])
                    if k in tot and tot[k] != v:
                        fail(f'{s} {plat} {r["type"]} {r["patient"]}: totals differ between pairs ({tot[k]} vs {v})')
                    tot[k] = v
# ---- Ti/Tv and het/hom per person: the four classes of one side summed
cls = {}
with open(R27, newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['stratum'] != 'genome':
            continue
        pair, plat = r['pair'].rsplit('_', 1)
        for s, sides in SIDES.items():
            if (pair, r['side']) == sides[0]:
                k = (s, plat, r['patient'])
                acc = cls.setdefault(k, {'n': 0, 'ti': 0, 'tv': 0, 'het': 0, 'hom': 0})
                for c in acc:
                    acc[c] += int(r[c])
rows = {}
for s in SETS:
    for plat in PLATFORMS:
        people = sorted({k[3] for k in tot if k[0] == s and k[1] == plat})
        if len(people) != 23:
            fail(f'{s} {plat}: {len(people)} individuals')
        for p in people:
            if (s, plat, p) not in cls:
                fail(f'{s} {plat}: no Ti/Tv row for one individual')
            if cls[(s, plat, p)]['n'] != tot[(s, plat, 'ALL', p)]:
                fail(f'{s} {plat}: the classes of r27_per_person.tsv do not add up to the total of all_counts.tsv for one individual')
        col = f'{plat.capitalize()} {s}'
        for typ, metric in (('ALL', 'Variant alleles per individual (M)'), ('SNV', 'SNVs (M)'), ('INDEL', 'Indels (M)')):
            rows.setdefault(metric, {})[col] = summary([tot[(s, plat, typ, p)] for p in people], dec=2, scale=1e6)
        rows.setdefault('Ti/Tv ratio (SNVs)', {})[col] = summary([cls[(s, plat, p)]['ti'] / cls[(s, plat, p)]['tv'] for p in people],
                                                                 ratio=True)
        rows.setdefault('Het/hom ratio', {})[col] = summary([cls[(s, plat, p)]['het'] / cls[(s, plat, p)]['hom'] for p in people],
                                                            ratio=True)

# ---- the FASTQ and BAM rows of the QC table (unchanged; DV used the GATK alignments)
out = []
with open(QC, newline='') as fh:
    qc = list(csv.reader(fh, delimiter='\t'))
head = qc[0]
for r in qc[1:]:
    if r[0].startswith('VCF-level'):
        break
    d = dict(zip(head, r + [''] * (len(head) - len(r))))
    if d['Illumina GATK'] == '' and d['Aviti GATK'] == '':
        out.append([r[0]] + [''] * len(COLS) + [''])          # a section header
        continue
    fastq = d['Illumina PPR'] == '' and d['Aviti PPR'] == ''
    vals = []
    for plat in ('Illumina', 'Aviti'):
        g, p = d[f'{plat} GATK'], d[f'{plat} PPR']
        vals += [g, '', ''] if fastq else [g, g, p]
    note = 'per platform: the same reads for the three pipelines' if fastq else 'DV used the GATK alignments'
    if r[0].startswith('Mean coverage'):
        note += '; PPR coverage includes duplicate reads (Supplementary Note 1)'
    out.append([r[0]] + [v.replace(chr(0x2013), '-') for v in vals] + [note])   # ranges with a hyphen
out.append(['VCF-level variant metrics (PASS)'] + [''] * len(COLS) + [''])
for metric in ('Variant alleles per individual (M)', 'SNVs (M)', 'Indels (M)', 'Ti/Tv ratio (SNVs)', 'Het/hom ratio'):
    out.append([metric] + [rows[metric][c] for c in COLS] + ['alleles as in the concordance analysis'])
path = os.path.join(HERE, 'Table1_revised.tsv')
with open(path, 'w', newline='', encoding='utf-8') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['Metric'] + COLS + ['Notes'])
    w.writerows(out)
with open(os.path.join(HERE, 'Table1_revised.provenance.txt'), 'w') as fh:
    for p in (QC, COUNTS, R27):
        fh.write(f'input\t{md5(p)}\t{p}\n')
    fh.write(f'script\t{md5(os.path.abspath(__file__))}\t{os.path.abspath(__file__)}\n')
    fh.write(f'table\t{md5(path)}\t{path}\n')
print(f'{len(out)} rows written to {path}')
print('TABLE1_REVISED_OK')
