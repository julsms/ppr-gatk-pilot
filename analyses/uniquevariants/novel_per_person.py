#!/usr/bin/env python3
"""novel_per_person.py - the novel records of one unique set in each individual (default GATK-only vs PPR; USET=PPRvsGATK, PPRvsDV or
DVvsPPR for the others). For each individual (numbered 1-23 in the order of the unique-variant figures), the number of unique records on
the autosomes, of novel ones, of novel SNVs and of novel T>G SNVs (T>G or A>C), Illumina and Aviti. Same novel definition as
05_tables.py: absent from dbSNP b157, gnomAD genomes and 1000 Genomes. per_person_summary.py gives the numbers of the Results: the share
of the three individuals with the most novel records and the range of the others.
Inputs : <UV_DIR>/carriers/<USET>_<platform>.multisample.vcf.gz (the unique records of each individual, from carriers.sh;
         read with bcftools query, only to count); <UV_DIR>/vep/vep_sites.tsv.gz (key, novel); default UV_DIR = uniquevariants
Outputs: per_person.tsv in the current folder (per_person_<USET>.tsv for the other sets; platform, individual, records, novel,
         novel_snv, novel_t_g, novel_het, novel_hom: genotype 0/1 or 1/0 = het, 1/1 or haploid 1 = hom, phased or not);
         per_person_ids.tsv (platform, individual, sample name); per_person.provenance.txt; last line NOVEL_PER_PERSON_DONE
Usage  : [USET=<set>] UV_DIR=<uniquevariants> python3 novel_per_person.py > novel_per_person.log
"""
import csv
import hashlib
import io
import os
import subprocess
import sys

import pandas as pd

U = os.environ.get('UV_DIR', 'uniquevariants')
USET = os.environ.get('USET', 'GATKvsPPR')
SUFFIX = '' if USET == 'GATKvsPPR' else f'_{USET}'
AUTO = {str(i) for i in range(1, 23)} | {f'chr{i}' for i in range(1, 23)}
T_G = {('T', 'G'), ('A', 'C')}


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


ids = []
inputs = [os.path.join(U, 'vep', 'vep_sites.tsv.gz')]
vep = pd.read_csv(inputs[0], sep='\t', dtype=str, keep_default_na=False, usecols=['key', 'novel'])
novel = set(vep.loc[vep['novel'] == '1', 'key'])
del vep
rows = []
for plat in ('illumina', 'aviti'):
    vcf = os.path.join(U, 'carriers', f'{USET}_{plat}.multisample.vcf.gz')
    inputs.append(vcf)
    names = subprocess.run(['bcftools', 'query', '-l', vcf], capture_output=True, text=True, check=True).stdout.split()
    if len(names) != 23:
        print(f'ERROR: {plat}: {len(names)} individuals, expected 23')
        print('NOVEL_PER_PERSON_FAILED')
        sys.exit(1)
    order = sorted(names, key=lambda c: (len(c), c))
    index = {c: order.index(c) + 1 for c in names}
    ids += [[plat, index[c], c] for c in order]
    out = subprocess.run(['bcftools', 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n', vcf], capture_output=True, text=True,
                         check=True).stdout
    q = pd.read_csv(io.StringIO(out), sep='\t', header=None, dtype=str, keep_default_na=False)
    q = q[q[0].isin(AUTO)]
    key = q[0] + ':' + q[1] + ':' + q[2] + ':' + q[3]
    is_novel = key.isin(novel)
    is_snv = (q[2].str.len() == 1) & (q[3].str.len() == 1)
    is_tg = pd.Series([(r, a) in T_G for r, a in zip(q[2], q[3])], index=q.index)
    for j, name in enumerate(names):
        gt = q[4 + j]
        has = ~gt.isin(['0/0', '0|0', './.', '.|.', '.', '0'])
        het = gt.isin(['0/1', '1/0', '0|1', '1|0'])
        hom = gt.isin(['1/1', '1|1', '1'])
        rows.append([plat, index[name], int(has.sum()), int((has & is_novel).sum()), int((has & is_novel & is_snv).sum()),
                     int((has & is_novel & is_snv & is_tg).sum()), int((het & is_novel).sum()), int((hom & is_novel).sum())])
    print(f'{plat}: {len(q)} autosomal records, {len(names)} individuals')
rows.sort(key=lambda r: (r[0] != 'illumina', r[1]))
with open(f'per_person{SUFFIX}.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['platform', 'individual', 'records', 'novel', 'novel_snv', 'novel_t_g', 'novel_het', 'novel_hom'])
    w.writerows(rows)
with open('per_person_ids.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['platform', 'individual', 'sample'])
    w.writerows(ids)
with open(f'per_person{SUFFIX}.provenance.txt', 'w') as fh:
    for p in inputs:
        fh.write(f'input\t{md5(p)}\t{p}\n')
    fh.write(f'script\t{md5(os.path.abspath(__file__))}\t{os.path.abspath(__file__)}\n')
print(f'{USET}: {len(rows)} rows written to per_person{SUFFIX}.tsv')
print('NOVEL_PER_PERSON_DONE')
