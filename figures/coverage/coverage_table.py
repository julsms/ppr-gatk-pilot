#!/usr/bin/env python3
"""coverage_table.py - data and correlations of Figure 6: per individual and platform the mean depth of the GATK alignments, the genotype
concordance of the call-set pairs and the k-mer support of the call sets; Spearman and Pearson correlations with depth, Holm-corrected
within each panel family (genome) and over the 10 regions. Inputs under PAPER_AGAIN_DIR (default .): qc/tables_GATKbam_qc/,
concordance/summary/all_counts.tsv, kmers_validation/revision/tables/per_sample_genome.tsv. Outputs next to this script.
"""
import csv
import hashlib
import os
import re
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
PA = os.environ.get('PAPER_AGAIN_DIR', '.')
QC = {p: os.path.join(PA, 'qc', 'tables_GATKbam_qc', f'bam_qc_per_sample_{p}_gatk.tsv') for p in ('illumina', 'aviti')}
COUNTS = os.path.join(PA, 'concordance', 'summary', 'all_counts.tsv')
MERFIN = os.path.join(PA, 'kmers_validation', 'revision', 'tables', 'per_sample_genome.tsv')
REGIONS = ['BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']
STRATA = ['genome'] + REGIONS
TYPES = ['ALL', 'SNV', 'INDEL']
PLATS = ['illumina', 'aviti']
WITHIN = ['GATK_DV', 'GATK_PPR', 'DV_PPR']
ACROSS = ['GATK', 'DV', 'PPR']
SETS = ['GATK_ALL', 'DV_ALL', 'PPR_ALL', 'GATKvsPPR', 'PPRvsGATK', 'DVvsPPR', 'PPRvsDV', 'PPRvsBOTH']


def md5(path):
    with open(path, 'rb') as fh:
        return hashlib.md5(fh.read()).hexdigest()


def fail(msg):
    print('ERROR:', msg)
    print('COVERAGE_TABLE_FAILED')
    sys.exit(1)


def read(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def patient_code(name):
    """the name of an individual in the concordance tables, from a QC sample name; numeric names pass unchanged"""
    m = re.match(r'^(\d+)ICDG(_S\d+)?$', name)
    return f'ICDG{m.group(1)}' if m else name


def holm(ps):
    order, adj, run = sorted(range(len(ps)), key=lambda i: ps[i]), [0.0] * len(ps), 0.0
    for k, i in enumerate(order):
        run = max(run, (len(ps) - k) * ps[i])
        adj[i] = min(1.0, run)
    return adj


# ---- concordance per patient (genotype concordance J_gt)
conc = {}
for r in read(COUNTS):
    if r['stratum'] in STRATA and r['type'] in TYPES:
        conc[(r['patient'], r['pair'], r['type'], r['stratum'])] = float(r['J_gt'])
codes = sorted({k[0] for k in conc})
if len(codes) != 23:
    fail(f'{len(codes)} individuals in all_counts.tsv, expected 23')
index = {c: i + 1 for i, c in enumerate(codes)}         # the numbering of concordance_table.py
pairs = [f'{p}_{pl}' for p in WITHIN for pl in PLATS] + [f'{s}_illumina_aviti' for s in ACROSS]
for c in codes:
    for p in pairs:
        for t in TYPES:
            for st in STRATA:
                if (c, p, t, st) not in conc:
                    fail(f'no concordance for one individual: {p} {t} {st}')
# ---- depth per patient and platform (GATK alignments)
depth = {}
for pl in PLATS:
    rows = read(QC[pl])
    if pl == 'aviti' and 'source_label' not in rows[0]:
        fail('the Aviti depth table has no source_label: not the table numbered by individual')
    flags = sorted({r['excluded'] for r in rows})
    print(f'{pl}: excluded values {flags}')
    for r in rows:
        if r['excluded'].strip().lower() in ('', 'no', 'false', '0'):
            depth[(patient_code(r['sample']), pl)] = float(r['mean_coverage'])
missing = [c for c in codes for pl in PLATS if (c, pl) not in depth]
if missing:
    fail(f'{len(missing)} individual-platform depths missing')
# ---- Merfin rate per patient, platform and call set (cohort, direct mode)
merf = {}
for r in read(MERFIN):
    if r['sample_kind'] == 'cohort' and r['mode'] == 'direct' and r['set_id'] in SETS:
        if r['sample'] not in index:
            fail('a Merfin sample is not one of the 23 individuals')
        merf[(r['sample'], r['platform'], r['set_id'])] = 100 * float(r['rate'])
done = {pl: sorted({k[0] for k in merf if k[1] == pl}) for pl in PLATS}
print('Merfin individuals done: ' + ', '.join(f'{pl} {len(done[pl])}' for pl in PLATS))

# ---- data table (numbered individuals)
data = []
for c in codes:
    i = index[c]
    for pl in PLATS:
        data.append([i, pl, 'depth', f'{depth[(c, pl)]:.3f}'])
        for s in SETS:
            if (c, pl, s) in merf:
                data.append([i, pl, f'merfin:{s}', f'{merf[(c, pl, s)]:.4f}'])
    for p in pairs:
        pl = 'both' if p.endswith('illumina_aviti') else p.rsplit('_', 1)[1]
        for t in TYPES:
            for st in STRATA:
                data.append([i, pl, f'concordance:{p.replace("_illumina_aviti", "_across") if pl == "both" else p.rsplit("_", 1)[0]}:{t}:{st}',
                             f'{100 * conc[(c, p, t, st)]:.4f}'])


# ---- correlations
def corr(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    rho, p = stats.spearmanr(x, y)
    r, pr = stats.pearsonr(x, y)
    return len(x), rho, p, r, pr


res = []                      # family, panel, x, y, platform, type, stratum, n, rho, p, r, pr
d_il = [depth[(c, 'illumina')] for c in codes]
d_av = [depth[(c, 'aviti')] for c in codes]
res.append(['depth', 'a', 'depth Illumina', 'depth Aviti', 'both', '', 'genome', *corr(d_il, d_av)])
for p in WITHIN:
    for pl in PLATS:
        for t in TYPES:
            for st in STRATA:
                res.append(['within', 'b' if (t, st) == ('ALL', 'genome') else 'supp', f'depth {pl}', p, pl, t, st,
                            *corr([depth[(c, pl)] for c in codes], [conc[(c, f'{p}_{pl}', t, st)] for c in codes])])
for s in ACROSS:
    for xname, xv in (('lower depth', [min(a, b) for a, b in zip(d_il, d_av)]), ('depth difference', [abs(a - b) for a, b in zip(d_il, d_av)])):
        for t in TYPES:
            for st in STRATA:
                res.append(['across', ('c' if xname == 'lower depth' else 'c2') if (t, st) == ('ALL', 'genome') else 'supp', xname,
                            f'{s} Illumina vs Aviti', 'both', t, st, *corr(xv, [conc[(c, f'{s}_illumina_aviti', t, st)] for c in codes])])
for s in SETS:
    for pl in PLATS:
        cs = done[pl]
        if len(cs) >= 5:
            res.append(['merfin', 'd', f'depth {pl}', s, pl, 'ALL', 'genome',
                        *corr([depth[(c, pl)] for c in cs], [merf[(c, pl, s)] for c in cs])])
# Holm: genome-level tests within each panel family; regional tests over the 10 regions of each (x, y, platform, type)
padj = [None] * len(res)
for panel in ('b', 'c', 'c2', 'd'):
    idx = [i for i, r in enumerate(res) if r[1] == panel]
    for i, a in zip(idx, holm([res[i][9] for i in idx])):
        padj[i] = a
groups = {}
for i, r in enumerate(res):
    if r[6] in REGIONS:
        groups.setdefault((r[0], r[2], r[3], r[4], r[5]), []).append(i)
for idx in groups.values():
    if len(idx) != len(REGIONS):
        fail('a regional Holm family without 10 regions')
    for i, a in zip(idx, holm([res[i][9] for i in idx])):
        padj[i] = a
for i, r in enumerate(res):
    if padj[i] is None:                 # the depth correlation and the genome rows of the supplement (SNV, INDEL): raw
        padj[i] = r[9]
with open(os.path.join(HERE, 'coverage_data.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['individual', 'platform', 'variable', 'value'])
    w.writerows(data)
with open(os.path.join(HERE, 'coverage_correlations.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['family', 'panel', 'x', 'y', 'platform', 'type', 'stratum', 'n', 'spearman_rho', 'p', 'p_holm', 'pearson_r', 'pearson_p'])
    for r, a in zip(res, padj):
        w.writerow(r[:8] + [f'{r[8]:.3f}', f'{r[9]:.3g}', f'{a:.3g}', f'{r[10]:.3f}', f'{r[11]:.3g}'])
for r, a in zip(res, padj):
    if r[1] in ('a', 'b', 'c', 'c2', 'd'):
        print(f'{r[1]:>2} {r[3]:<26} {r[4]:<9} n {r[7]:>2}  rho {r[8]:+.3f} p {r[9]:.2g} holm {a:.2g}  r {r[10]:+.3f}')
with open(os.path.join(HERE, 'coverage_table.provenance.txt'), 'w') as fh:
    for p in (QC['illumina'], QC['aviti'], COUNTS, MERFIN):
        fh.write(f'input\t{md5(p)}\t{p}\n')
    fh.write(f'script\t{md5(os.path.abspath(__file__))}\t{os.path.abspath(__file__)}\n')
    for name in ('coverage_data.tsv', 'coverage_correlations.tsv'):
        fh.write(f'table\t{md5(os.path.join(HERE, name))}\t{os.path.join(HERE, name)}\n')
print(f'{len(data)} data rows, {len(res)} correlations')
print('COVERAGE_TABLE_OK')
