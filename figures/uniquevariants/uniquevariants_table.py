#!/usr/bin/env python3
"""uniquevariants_table.py - the data of Figure 4 and Supplementary Figures S4 and S13-S18, from the outputs of analyses/uniquevariants
(UV_DIR): unique variants per individual (numbered 1-23), frequencies in gnomAD genomes and 1000 Genomes (largest value of a field,
present = AF > 0), consequences, GATK status, Ti/Tv and het/hom, singletons (RELATED left out), HG002 precision, novel alleles.
alleles = distinct alleles of the cohort; records = alleles counted once per individual. Outputs: uv_*.tsv in this folder.
"""
import csv
import hashlib
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
U = os.environ.get('UV_DIR', 'uniquevariants')
# one of the two first-degree relatives, left out of the singletons
RELATED = os.environ.get('RELATED', '6')
SETS =['GATKvsPPR', 'PPRvsGATK', 'DVvsPPR', 'PPRvsDV', 'PPRvsBOTH']
PLATFORMS = ['illumina', 'aviti']
REGIONS = ['BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']
STRATA = ['genome'] + REGIONS
TYPES = ['ALL', 'SNV', 'INDEL']
SOURCES = [('gnomADg', 'gnomADg_AF'), ('1000G', 'AF')]
EDGES = np.round(np.arange(-6.0, 0.0001, 0.25), 2)
CODING = [('missense', ['missense_variant', 'protein_altering_variant']),
          ('synonymous', ['synonymous_variant', 'stop_retained_variant']),
          ('splice', ['splice_acceptor_variant', 'splice_donor_variant', 'splice_region_variant', 'splice_donor_5th_base_variant',
                      'splice_donor_region_variant', 'splice_polypyrimidine_tract_variant']),
          ('frameshift', ['frameshift_variant']), ('in-frame', ['inframe_insertion', 'inframe_deletion']),
          ('stop gained, stop lost, start lost', ['stop_gained', 'stop_lost', 'start_lost'])]
CODING_ALL = [m for _, ms in CODING for m in ms] + ['coding_sequence_variant', 'incomplete_terminal_codon_variant',
                                                     'start_retained_variant']
# five groups that add up to every allele: the most severe consequence of the allele (consequences.tsv) in one of them
COMPOSITION = [('intergenic', ['intergenic_variant']), ('intron', ['intron_variant']),
               ('upstream or downstream of a gene', ['upstream_gene_variant', 'downstream_gene_variant']),
               ('coding or splice', CODING_ALL)]
OTHER_GROUP = 'UTR, non-coding exon, regulatory or other'
AF_CLASSES = [('lt0.1pct', 0.0, 0.001), ('0.1-1pct', 0.001, 0.01), ('1-5pct', 0.01, 0.05), ('ge5pct', 0.05, 1.01)]
R119 = [('vqsr', ['a_vqsr']), ('excesshet', ['a2_excesshet']), ('other', ['a3_other_filter', 'b_absent', 'a_unmatched']),
        ('ref_in_person', ['c_ref']), ('missing_in_person', ['c_missing']), ('not_called', ['d_none'])]


def md5(path):
    with open(path, 'rb') as fh:
        return hashlib.md5(fh.read()).hexdigest()


def fail(msg):
    print('ERROR:', msg)
    print('UV_TABLE_FAILED')
    sys.exit(1)


def max_af(v):
    best = np.nan
    for x in str(v).split('&'):
        try:
            f = float(x)
        except ValueError:
            continue
        if f > 0 and not f <= best:
            best = f
    return best


def write(name, head, rows):
    path = os.path.join(HERE, name)
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(head)
        w.writerows(rows)
    print(f'{len(rows)} rows written to {path}')
    return path


inputs = []

# ---- unique variants per individual
src = os.path.join(U, 'tables', 'unique_counts.tsv')
inputs.append(src)
uc = pd.read_csv(src, sep='\t', dtype=str, keep_default_na=False)
uc = uc[uc['stratum'].isin(STRATA) & uc['type'].isin(TYPES)]
uc[['set', 'platform']] = uc['uset'].str.rsplit('_', n=1, expand=True)
uc = uc[uc['set'].isin(SETS)]
codes = sorted(uc['patient'].unique(), key=lambda c: (len(c), c))
if len(codes) != 23:
    fail(f'{len(codes)} individuals, expected 23')
index = {c: i + 1 for i, c in enumerate(codes)}
for (s, p, st, t), g in uc.groupby(['set', 'platform', 'stratum', 'type']):
    if g['patient'].nunique() != 23 or len(g) != 23:
        fail(f'{s} {p} {st} {t}: {len(g)} rows, expected one per individual')
if len(uc) != 23 * len(SETS) * len(PLATFORMS) * len(STRATA) * len(TYPES):
    fail(f'{len(uc)} count rows, expected {23 * len(SETS) * len(PLATFORMS) * len(STRATA) * len(TYPES)}')
counts = [[index[r.patient], r.set, r.platform, r.stratum, r.type, int(r.n)] for r in uc.itertuples()]
counts.sort(key=lambda x: (SETS.index(x[1]), PLATFORMS.index(x[2]), STRATA.index(x[3]), TYPES.index(x[4]), x[0]))
genome_records = {(r[1], r[2]): 0 for r in counts}
for r in counts:
    if r[3] == 'genome' and r[4] == 'ALL':
        genome_records[(r[1], r[2])] += r[5]
if RELATED not in index:
    fail(f'RELATED {RELATED} is not one of the individuals')
genome_records_unrel = {k: 0 for k in genome_records}
for r in counts:
    if r[3] == 'genome' and r[4] == 'ALL' and r[0] != index[RELATED]:
        genome_records_unrel[(r[1], r[2])] += r[5]

# ---- population frequencies of the unique alleles
vep_path = os.path.join(U, 'vep', 'vep_sites.tsv.gz')
inputs.append(vep_path)
vep = pd.read_csv(vep_path, sep='\t', dtype=str, keep_default_na=False, usecols=['key'] + [c for _, c in SOURCES])
for _, c in SOURCES:
    vep[c] = vep[c].map(max_af)
vep = vep.drop_duplicates('key').set_index('key')
hist, af_classes = [], []
carriers = {}
for s in SETS:
    for p in PLATFORMS:
        cp = os.path.join(U, 'carriers', f'{s}_{p}.carriers.tsv.gz')
        inputs.append(cp)
        c = pd.read_csv(cp, sep='\t', dtype=str, keep_default_na=False, usecols=['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'n_U'] + REGIONS)
        c['key'] = c['CHROM'] + ':' + c['POS'] + ':' + c['REF'] + ':' + c['ALT']
        c['n_U'] = c['n_U'].astype(int)
        miss = int((~c['key'].isin(vep.index)).sum())
        if miss:
            fail(f'{s} {p}: {miss} alleles without a VEP row')
        if int(c['n_U'].sum()) != genome_records[(s, p)]:
            fail(f'{s} {p}: {int(c["n_U"].sum())} records in the carrier table, {genome_records[(s, p)]} in unique_counts.tsv')
        j = c.join(vep, on='key')
        carriers[(s, p)] = len(c)
        for st in STRATA:
            ms = np.ones(len(j), bool) if st == 'genome' else (j[st] == '1').to_numpy()
            for t in TYPES:
                m = ms if t == 'ALL' else ms & (j['TYPE'] == t).to_numpy()
                sub = j[m]
                rec = sub['n_U'].to_numpy()
                for src_name, col in SOURCES:
                    af = sub[col].to_numpy(dtype=float)
                    present = np.isfinite(af)
                    la = np.log10(af[present])
                    if len(la) and (la.min() < EDGES[0] or la.max() > 0):
                        fail(f'{s} {p} {st} {t} {src_name}: log10 AF outside {EDGES[0]} to 0')
                    nk, _ = np.histogram(la, bins=EDGES)
                    nr, _ = np.histogram(la, bins=EDGES, weights=rec[present])
                    if int((~present).sum() + nk.sum()) != len(sub):
                        fail(f'{s} {p} {st} {t} {src_name}: the bins do not add up to the alleles')
                    hist.append([s, p, st, t, src_name, 'absent', 'absent', int((~present).sum()), int(rec[~present].sum())])
                    for lo, hi, a, b in zip(EDGES[:-1], EDGES[1:], nk, nr):
                        hist.append([s, p, st, t, src_name, f'{lo:.2f}', f'{hi:.2f}', int(a), int(b)])
                    af_classes.append([s, p, st, t, src_name, 'absent', int((~present).sum()), int(rec[~present].sum())])
                    in_class = 0
                    for name, lo, hi in AF_CLASSES:
                        mk = present.copy()
                        mk[present] = (af[present] >= lo) & (af[present] < hi)
                        in_class += int(mk.sum())
                        af_classes.append([s, p, st, t, src_name, name, int(mk.sum()), int(rec[mk].sum())])
                    if in_class != int(present.sum()):
                        fail(f'{s} {p} {st} {t} {src_name}: the frequency classes do not add up to the alleles found')
        print(f'{s} {p}: {len(c)} alleles, {int(c["n_U"].sum())} records')

# ---- coding and splice consequences (distinct alleles of the cohort)
src = os.path.join(U, 'tables', 'consequences.tsv')
inputs.append(src)
cq = pd.read_csv(src, sep='\t', dtype=str, keep_default_na=False)
cq[['set', 'platform']] = cq['uset'].str.rsplit('_', n=1, expand=True)
cq = cq[cq['set'].isin(SETS) & cq['stratum'].isin(STRATA)]
cq['n'] = cq['n'].astype(int)
cons = []
for (s, p, st), g in cq.groupby(['set', 'platform', 'stratum']):
    n_region = set(g['n_stratum'])
    if len(n_region) != 1:
        fail(f'{s} {p} {st}: more than one region total')
    n_region = int(n_region.pop())
    if g['n'].sum() != n_region:
        fail(f'{s} {p} {st}: the consequence classes add up to {g["n"].sum()}, not {n_region}')
    if st == 'genome' and n_region != carriers[(s, p)]:
        fail(f'{s} {p}: {n_region} alleles in consequences.tsv, {carriers[(s, p)]} in the carrier table')
    for name, members in CODING:
        cons.append([s, p, st, 'coding', name, int(g.loc[g['consequence_class'].isin(members), 'n'].sum()), n_region])
    used = 0
    for name, members in COMPOSITION:
        v = int(g.loc[g['consequence_class'].isin(members), 'n'].sum())
        used += v
        cons.append([s, p, st, 'composition', name, v, n_region])
    cons.append([s, p, st, 'composition', OTHER_GROUP, n_region - used, n_region])
GROUP_ORDER = [c for c, _ in CODING] + [c for c, _ in COMPOSITION] + [OTHER_GROUP]
cons.sort(key=lambda x: (SETS.index(x[0]), PLATFORMS.index(x[1]), STRATA.index(x[2]), x[3], GROUP_ORDER.index(x[4])))
n_cons = len(SETS) * len(PLATFORMS) * len(STRATA) * (len(CODING) + len(COMPOSITION) + 1)
if len(cons) != n_cons:
    fail(f'{len(cons)} consequence rows, expected {n_cons}')

# ---- why GATK did not report the PPR-only variants (classes against the raw GATK calls)
src = os.path.join(U, 'tables', 'r119_summary.tsv')
inputs.append(src)
rs = pd.read_csv(src, sep='\t', dtype=str, keep_default_na=False)
rs = rs[(rs['mode'] == 'new') & rs['stratum'].isin(STRATA) & rs['type'].isin(TYPES)]
r119 = []
for r in rs.itertuples():
    n = int(r.n)
    parts = [(name, sum(int(getattr(r, k)) for k in keys)) for name, keys in R119]
    if sum(v for _, v in parts) != n:
        fail(f'GATK status {r.platform} {r.stratum} {r.type}: the classes add up to {sum(v for _, v in parts)}, not {n}')
    if r.stratum == 'genome' and r.type == 'ALL' and n != genome_records[('PPRvsGATK', r.platform)]:
        fail(f'GATK status {r.platform}: {n} PPR-only records, {genome_records[("PPRvsGATK", r.platform)]} in unique_counts.tsv')
    for name, v in parts:
        r119.append([r.platform, r.stratum, r.type, name, v, f'{v / n:.6f}' if n else 'NA'])
if len(r119) != len(PLATFORMS) * len(STRATA) * len(TYPES) * len(R119):
    fail(f'{len(r119)} GATK status rows, expected {len(PLATFORMS) * len(STRATA) * len(TYPES) * len(R119)}')

# ---- Ti/Tv and het/hom of the shared and unique variants (cohort sums)
src = os.path.join(U, 'tables', 'r27_titv_hethom.tsv')
inputs.append(src)
tt = pd.read_csv(src, sep='\t', dtype=str, keep_default_na=False)
tt = tt[tt['class'].isin(['shared', 'unique']) & tt['stratum'].isin(STRATA)]
titv = []
for r in tt.to_dict('records'):
    snv, ti, tv, het, hom = (int(r['snv']), int(r['ti']), int(r['tv']), int(r['het']), int(r['hom']))
    where = f'Ti/Tv and het/hom {r["pair"]} {r["side"]} {r["class"]} {r["stratum"]}'
    if ti + tv != snv:
        fail(f'{where}: Ti + Tv is not the SNV count')
    for mine, theirs, what in ((ti / tv if tv else float('nan'), r['titv'], 'Ti/Tv'),
                               (het / hom if hom else float('nan'), r['hethom'], 'het/hom')):
        if theirs not in ('', 'NA', 'nan') and abs(mine - float(theirs)) > 1e-4:
            fail(f'{where}: {what} {mine:.6f} here, {theirs} in the table')
    pair, plat = r['pair'].rsplit('_', 1)
    titv.append([pair, plat, r['side'], r['set'], r['class'], r['stratum'], snv, ti, tv, het, hom, f'{ti / tv:.4f}' if tv else 'NA',
                 f'{het / hom:.4f}' if hom else 'NA'])
if len(titv) != 3 * len(PLATFORMS) * 2 * 2 * len(STRATA):
    fail(f'{len(titv)} Ti/Tv and het/hom rows, expected {3 * len(PLATFORMS) * 2 * 2 * len(STRATA)} (3 pairs, 2 sides, shared and unique)')

# ---- singletons: the allele is carried by no other individual in that pipeline's calls, among the 22 individuals other than
# RELATED (05_tables.py: n_excl_related and private_full_excl_related); RELATED counts as not sequenced: its own unique records are
# left out, and its calls do not make the allele of another individual shared
src = os.path.join(U, 'tables', 'singletons_summary.tsv')
inputs.append(src)
sg = pd.read_csv(src, sep='\t', dtype=str, keep_default_na=False)
sg[['set', 'platform']] = sg['uset'].str.rsplit('_', n=1, expand=True)
sg = sg[sg['set'].isin(SETS) & sg['stratum'].isin(STRATA) & sg['type'].isin(TYPES)]
sing = []
for r in sg.itertuples():
    n, k = int(r.n_excl_related), int(r.private_full_excl_related)
    if (r.stratum, r.type) == ('genome', 'ALL') and n != genome_records_unrel[(r.set, r.platform)]:
        fail(f'singletons {r.set} {r.platform}: {n} records, {genome_records_unrel[(r.set, r.platform)]} in unique_counts.tsv '
             f'without {RELATED}')
    sing.append([r.set, r.platform, r.stratum, r.type, n, k, f'{k / n:.6f}' if n else 'NA'])
if len(sing) != len(SETS) * len(PLATFORMS) * len(STRATA) * len(TYPES):
    fail(f'{len(sing)} singleton rows, expected {len(SETS) * len(PLATFORMS) * len(STRATA) * len(TYPES)}')

# ---- HG002: precision of the unique calls against GIAB v5.0q, inside BB (the genome row holds the calls inside BB)
prec = []
for p in PLATFORMS:
    src = os.path.join(U, 'hg002', 'tables', f'precision_{p}.tsv')
    inputs.append(src)
    t = pd.read_csv(src, sep='\t', dtype=str, keep_default_na=False)
    t = t[(t['class'] == 'unique') & (t['consequence'] == 'ALL') & t['stratum'].isin(STRATA) & t['type'].isin(TYPES)]
    for r in t.itertuples():
        s = r.uset.rsplit('_', 1)[0]
        if s not in SETS:
            continue
        tp, fp, fg, ne, n = int(r.TP), int(r.FP), int(r.FP_gt), int(r.NE), int(r.n)
        if tp + fp + fg + ne != n:
            fail(f'HG002 {s} {p} {r.stratum} {r.type}: TP + FP + FP_gt + outside BB is not n')
        pr = tp / (tp + fp + fg) if tp + fp + fg else float('nan')
        if r.precision not in ('', 'NA') and abs(pr - float(r.precision)) > 1e-6:
            fail(f'HG002 {s} {p} {r.stratum} {r.type}: precision {pr:.6f} here, {r.precision} in the table')
        prec.append([s, p, r.stratum, r.type, n, tp, fp, fg, ne, f'{pr:.6f}' if tp + fp + fg else 'NA'])
if len(prec) != len(SETS) * len(PLATFORMS) * len(STRATA) * len(TYPES):
    fail(f'{len(prec)} HG002 rows, expected {len(SETS) * len(PLATFORMS) * len(STRATA) * len(TYPES)}')

# ---- novel alleles: absent from dbSNP b157, gnomAD genomes and 1000 Genomes; distinct alleles of the cohort on the autosomes
# (05_tables.py, novel.tsv)
src = os.path.join(U, 'tables', 'novel.tsv')
inputs.append(src)
nv = pd.read_csv(src, sep='\t', dtype=str, keep_default_na=False)
nv[['set', 'platform']] = nv['uset'].str.rsplit('_', n=1, expand=True)
NOVEL_STRATA = ['autosomes'] + REGIONS
nv = nv[nv['set'].isin(SETS) & nv['stratum'].isin(NOVEL_STRATA)]
novel = []
for r in nv.itertuples():
    n, k = int(r.n), int(r.novel)
    if k > n:
        fail(f'novel {r.uset} {r.stratum}: {k} novel alleles of {n}')
    if n and r.share_novel not in ('', 'NA') and abs(k / n - float(r.share_novel)) > 1e-6:
        fail(f'novel {r.uset} {r.stratum}: share {k / n:.6f} here, {r.share_novel} in the table')
    novel.append([r.set, r.platform, r.stratum, n, k, f'{k / n:.6f}' if n else 'NA'])
if len(novel) != len(SETS) * len(PLATFORMS) * len(NOVEL_STRATA):
    fail(f'{len(novel)} novel rows, expected {len(SETS) * len(PLATFORMS) * len(NOVEL_STRATA)}')

outs = [write('uv_counts.tsv', ['individual', 'set', 'platform', 'region', 'type', 'n'], counts),
        write('uv_af_hist.tsv', ['set', 'platform', 'region', 'type', 'source', 'bin_lo', 'bin_hi', 'alleles', 'records'], hist),
        write('uv_af_classes.tsv', ['set', 'platform', 'region', 'type', 'source', 'class', 'alleles', 'records'], af_classes),
        write('uv_consequences.tsv', ['set', 'platform', 'region', 'kind', 'group', 'alleles', 'alleles_in_region'], cons),
        write('uv_r119.tsv', ['platform', 'region', 'type', 'class', 'records', 'share'], r119),
        write('uv_titv_hethom.tsv', ['pair', 'platform', 'side', 'set', 'class', 'region', 'snv', 'ti', 'tv', 'het', 'hom', 'titv',
                                     'hethom'], titv),
        write('uv_singletons.tsv', ['set', 'platform', 'region', 'type', 'records', 'singletons', 'share'], sing),
        write('uv_hg002_precision.tsv', ['set', 'platform', 'region', 'type', 'n', 'TP', 'FP', 'FP_gt', 'outside_BB', 'precision'],
              prec),
        write('uv_novel.tsv', ['set', 'platform', 'region', 'alleles', 'novel', 'share'], novel)]
with open(os.path.join(HERE, 'uv_data.provenance.txt'), 'w') as fh:
    for path in inputs:
        fh.write(f'input\t{md5(path)}\t{path}\n')
    fh.write(f'script\t{md5(os.path.abspath(__file__))}\t{os.path.abspath(__file__)}\n')
    for path in outs:
        fh.write(f'table\t{md5(path)}\t{path}\n')
print('UV_TABLE_OK')
