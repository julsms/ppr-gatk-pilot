#!/usr/bin/env python3
"""novel_profile.py - the novel alleles of each unique set on the autosomes (05_tables.py: novel = absent from dbSNP b157, gnomAD
genomes and 1000 Genomes; distinct alleles of the cohort), by variant type, by substitution class of the SNVs, by the number of
individuals who carry the allele in the full call set of the pipeline, by region, and by the database each allele is absent from.
The Results quote it for the GATK-only alleles on Aviti (share of SNVs, T>G or A>C, one carrier).
Inputs : <UV_DIR>/carriers/<set>_<platform>.carriers.tsv.gz (distinct alleles of each unique set: TYPE, AUTO, region flags,
         n_carriers_full = individuals with the allele in the full call set, n_U = individuals where it is unique);
         <UV_DIR>/vep/vep_sites.tsv.gz (key, in_dbSNP, in_gnomADg, in_1000G, novel); default UV_DIR = uniquevariants.
         Read only; no VCF is read.
Outputs: in the current folder: profile_overview.tsv, profile_types.tsv, profile_substitutions.tsv, profile_carriers.tsv,
         profile_regions.tsv, profile.provenance.txt (md5 of the inputs and of this script); last line NOVEL_PROFILE_DONE
Usage  : UV_DIR=<uniquevariants> python3 novel_profile.py > novel_profile.log
"""
import csv
import hashlib
import os
import sys

import pandas as pd

U = os.environ.get('UV_DIR', 'uniquevariants')
USETS = [f'{s}_{p}' for s in ('GATKvsPPR', 'PPRvsGATK', 'DVvsPPR', 'PPRvsDV') for p in ('illumina', 'aviti')]
REGIONS = ['BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']
COMP = {('A', 'C'): 'T>G', ('A', 'G'): 'T>C', ('A', 'T'): 'T>A', ('C', 'A'): 'C>A', ('C', 'G'): 'C>G', ('C', 'T'): 'C>T',
        ('G', 'A'): 'C>T', ('G', 'C'): 'C>G', ('G', 'T'): 'C>A', ('T', 'A'): 'T>A', ('T', 'C'): 'T>C', ('T', 'G'): 'T>G'}
BINS = [(1, 1, '1'), (2, 2, '2'), (3, 5, '3-5'), (6, 10, '6-10'), (11, 23, '11-23')]


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write(name, head, rows):
    with open(name, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(head)
        w.writerows(rows)
    print(f'{len(rows)} rows written to {name}')


def pct(a, b):
    return f'{100 * a / b:.2f}' if b else 'NA'


inputs = [os.path.join(U, 'vep', 'vep_sites.tsv.gz')] + [os.path.join(U, 'carriers', f'{u}.carriers.tsv.gz') for u in USETS]
for p in inputs:
    if not os.path.isfile(p):
        print('ERROR: missing input', p)
        print('NOVEL_PROFILE_FAILED')
        sys.exit(1)
vep = pd.read_csv(inputs[0], sep='\t', dtype=str, keep_default_na=False,
                  usecols=['key', 'in_dbSNP', 'in_gnomADg', 'in_1000G', 'novel'])
for c in ('in_dbSNP', 'in_gnomADg', 'in_1000G', 'novel'):
    vep[c] = vep[c] == '1'
print(f'{len(vep)} VEP sites')
overview, types, subs, carriers, regions = [], [], [], [], []
for u, path in zip(USETS, inputs[1:]):
    cols = ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'AUTO', 'n_carriers_full', 'n_U'] + REGIONS
    k = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, usecols=cols)
    k = k[k['AUTO'] == '1'].copy()
    k['key'] = k['CHROM'] + ':' + k['POS'] + ':' + k['REF'] + ':' + k['ALT']
    j = k.merge(vep, on='key', how='left')
    if j['novel'].isna().any():
        print(f'ERROR: {u}: {int(j["novel"].isna().sum())} alleles without a VEP row')
        print('NOVEL_PROFILE_FAILED')
        sys.exit(1)
    j['ncar'] = j['n_carriers_full'].astype(int)
    nov = j[j['novel']]
    n, nn = len(j), len(nov)
    overview.append([u, n, nn, pct(nn, n), int((~j['in_dbSNP']).sum()), int((~j['in_gnomADg']).sum()), int((~j['in_1000G']).sum()),
                     int((nov['BB'] == '1').sum()), pct(int((nov['BB'] == '1').sum()), nn), int((nov['ncar'] == 1).sum()),
                     pct(int((nov['ncar'] == 1).sum()), nn), int((j['ncar'] == 1).sum()), pct(int((j['ncar'] == 1).sum()), n)])
    for t, g in j.groupby('TYPE'):
        gn = int(g['novel'].sum())
        types.append([u, t, len(g), gn, pct(gn, len(g)), pct(gn, nn)])
    snv = j[j['TYPE'] == 'SNV'].copy()
    snv['sub'] = [COMP.get((r, a), 'other') for r, a in zip(snv['REF'], snv['ALT'])]
    for s, g in snv.groupby('sub'):
        gn = int(g['novel'].sum())
        subs.append([u, s, len(g), pct(len(g), len(snv)), gn, pct(gn, int(snv['novel'].sum()))])
    for lo, hi, lab in BINS:
        g = j[(j['ncar'] >= lo) & (j['ncar'] <= hi)]
        gn = int(g['novel'].sum())
        carriers.append([u, lab, len(g), gn, pct(gn, len(g)), pct(gn, nn)])
    for r in REGIONS:
        g = j[j[r] == '1']
        gn = int(g['novel'].sum())
        regions.append([u, r, len(g), gn, pct(gn, len(g))])
    print(f'{u}: {n} autosomal alleles, {nn} novel')
write('profile_overview.tsv', ['uset', 'alleles', 'novel', 'novel_pct', 'absent_dbSNP', 'absent_gnomADg', 'absent_1000G',
                               'novel_in_BB', 'novel_in_BB_pct', 'novel_one_carrier', 'novel_one_carrier_pct', 'alleles_one_carrier',
                               'alleles_one_carrier_pct'], overview)
write('profile_types.tsv', ['uset', 'type', 'alleles', 'novel', 'novel_pct_of_type', 'pct_of_novel'], types)
write('profile_substitutions.tsv', ['uset', 'substitution', 'snvs', 'pct_of_snvs', 'novel_snvs', 'pct_of_novel_snvs'], subs)
write('profile_carriers.tsv', ['uset', 'carriers_full', 'alleles', 'novel', 'novel_pct_of_bin', 'pct_of_novel'], carriers)
write('profile_regions.tsv', ['uset', 'region', 'alleles', 'novel', 'novel_pct'], regions)
with open('profile.provenance.txt', 'w') as fh:
    for p in inputs:
        fh.write(f'input\t{md5(p)}\t{p}\n')
    fh.write(f'script\t{md5(os.path.abspath(__file__))}\t{os.path.abspath(__file__)}\n')
print('NOVEL_PROFILE_DONE')
