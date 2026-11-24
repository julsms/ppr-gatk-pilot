#!/usr/bin/env python3
"""Data of Figure 7: recall, precision, F1, FN and FP of every GATK call set on HG002 (GIAB v5.0q) and HG001 (GIAB v4.2.1), by
region, type and filter, with SNVs and indels together as type ALL; the comparisons quoted in the text on stdout.
Inputs: data/HG002/<set>.summary.tsv and data/HG001/<set>.summary.tsv next to this script (hap.py summaries, see README.md).
Output: figure7_table.tsv. Usage: python3 figure7_table.py > figure7_table.txt
"""
import csv
import glob
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'figure7_table.tsv')


def label(genome, name):
    """(platform, calling, filter label) of a set folder name"""
    if genome == 'HG002':
        p, rest = name.split('_', 1)
        rest = rest.replace('GATK_HaplotypeCaller_', '')
        if rest in ('joint', 'single'):
            return p, 'cohort' if rest == 'joint' else 'single', 'pipeline'
        m = re.match(r'joint_snp([\d.]+)_indel([\d.]+)$', rest)
        if m:
            return p, 'cohort', f'pipeline levels {m.group(1)}/{m.group(2)}'
        m = re.match(r'joint_vqsrbp_snp([\d.]+)_indel([\d.]+)$', rest)
        if m:
            return p, 'cohort', f'BP {m.group(1)}/{m.group(2)}'
        return p, 'DeepVariant', rest
    m = re.match(r'HG001_(\w+?)_GATK_HaplotypeCaller_(cohort|single)_(.*)$', name)
    p, calling, rest = m.group(1), m.group(2), m.group(3)
    if rest == 'pipeline_vqsr':
        return p, calling, 'pipeline'
    m2 = re.match(r'vqsrbp(_no_excesshet)?_snp([\d.]+)_indel([\d.]+)$', rest)
    return p, calling, f'BP {m2.group(2)}/{m2.group(3)}' + (' no ExcessHet' if m2.group(1) else '')


rows = []
for genome in ('HG002', 'HG001'):
    for f in sorted(glob.glob(os.path.join(HERE, 'data', genome, '*.summary.tsv'))):
        name = os.path.basename(f)[:-len('.summary.tsv')]
        p, calling, filt = label(genome, name)
        both = {}
        with open(f, newline='') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                fn, fp = int(r['FN']), int(r['FP'])
                rows.append({'genome': genome, 'platform': p, 'calling': calling, 'setting': filt, 'set': name, 'region': r['region'],
                             'type': r['type'], 'filter': r['filter'], 'recall': float(r['recall']), 'precision': float(r['precision']),
                             'F1': float(r['F1']), 'FN': fn, 'FP': fp, 'errors': fn + fp})
                both.setdefault((r['region'], r['filter']), []).append(r)
        # ALL = SNVs and indels together: recall = (TP + TP) / (truth + truth); precision from the query true
        # positives, which hap.py does not print in summary.tsv but which are query_total - FP - query_unk (checked: this gives back
        # the printed precision of each type)
        for (region, filt_), pair in both.items():
            if len(pair) != 2:
                continue
            tp = sum(int(r['TP']) for r in pair)
            truth = sum(int(r['truth_total']) for r in pair)
            fp = sum(int(r['FP']) for r in pair)
            fn = sum(int(r['FN']) for r in pair)
            qtp = sum(int(r['query_total']) - int(r['FP']) - int(r['query_unk']) for r in pair)
            rec, prec = tp / truth, qtp / (qtp + fp)
            rows.append({'genome': genome, 'platform': p, 'calling': calling, 'setting': filt, 'set': name, 'region': region,
                         'type': 'ALL', 'filter': filt_, 'recall': rec, 'precision': prec, 'F1': 2 * rec * prec / (rec + prec),
                         'FN': fn, 'FP': fp, 'errors': fn + fp})
# check of the query true positives: the precision computed from them equals the printed one for every SNP and INDEL row
worst = max(abs((int(r['query_total']) - int(r['FP']) - int(r['query_unk'])) / (int(r['query_total']) - int(r['query_unk'])) - float(r['precision']))
            for f in glob.glob(os.path.join(HERE, 'data', '*', '*.summary.tsv')) for r in csv.DictReader(open(f), delimiter='\t')
            if int(r['query_total']) - int(r['query_unk']) > 0)
with open(OUT, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(rows)


def get(genome, p, calling, setting, region='BB', t='SNP', filt='PASS'):
    for r in rows:
        if (r['genome'], r['platform'], r['calling'], r['setting'], r['region'], r['type'], r['filter']) == (genome, p, calling, setting, region, t, filt):
            return r
    return None


def fmt(r):
    return f"F1 {r['F1']:.4f} recall {r['recall']:.4f} precision {r['precision']:.4f} errors {r['errors']:,} (FN {r['FN']:,} FP {r['FP']:,})"


print(f'== check: largest difference between the precision from query_total - FP - query_unk and the printed precision: {worst:.2e}')
print('== HG001: cohort vs alone, all variants (SNVs and indels together), PASS, by region, at each Best Practices level pair')
for p in sorted({r['platform'] for r in rows if r['genome'] == 'HG001'}):
    for s in ('BP 99.0/99.0', 'BP 99.9/99.0', 'BP 99.9/99.9'):
        for region in ('BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG'):
            a, b = get('HG001', p, 'single', s, region, 'ALL'), get('HG001', p, 'cohort', s, region, 'ALL')
            if a and b:
                print(f'HG001 {p} [{s}] {region} ALL: alone F1 {a["F1"]:.4f} errors {a["errors"]:,}; cohort F1 {b["F1"]:.4f} errors '
                      f'{b["errors"]:,} ({100 * (b["errors"] - a["errors"]) / a["errors"]:+.1f}%)')

print('== every set, PASS, BB and LM')
seen = sorted({(r['genome'], r['platform'], r['calling'], r['setting']) for r in rows})
for g, p, c, s in seen:
    for region in ('BB', 'LM'):
        for t in ('SNP', 'INDEL'):
            r = get(g, p, c, s, region, t)
            if r:
                print(f'{g} {p} {c} [{s}] {region} {t}: {fmt(r)}')

print('== INDEL 99.9 vs 99.0 at SNP 99.9 (Best Practices): change of errors, PASS (negative = fewer with 99.9)')
for g, p, c in sorted({(g, p, c) for g, p, c, s in seen if s == 'BP 99.9/99.9'}):
    for region in ('BB', 'LM', 'TRHP', 'AD'):
        for t in ('SNP', 'INDEL'):
            a, b = get(g, p, c, 'BP 99.9/99.0', region, t), get(g, p, c, 'BP 99.9/99.9', region, t)
            if a and b:
                print(f'{g} {p} {c} {region} {t}: errors {a["errors"]:,} -> {b["errors"]:,} ({b["errors"] - a["errors"]:+,}; '
                      f'{100 * (b["errors"] - a["errors"]) / a["errors"]:+.1f}%), F1 {a["F1"]:.4f} -> {b["F1"]:.4f}')

print('== SNP 99.9 vs 99.0 at INDEL 99.0 (Best Practices; HG001 only has both)')
for g, p, c in sorted({(g, p, c) for g, p, c, s in seen if s == 'BP 99.0/99.0'}):
    for region in ('BB', 'LM'):
        a, b = get(g, p, c, 'BP 99.0/99.0', region, 'SNP'), get(g, p, c, 'BP 99.9/99.0', region, 'SNP')
        if a and b:
            print(f'{g} {p} {c} {region} SNP: errors {a["errors"]:,} -> {b["errors"]:,} ({100 * (b["errors"] - a["errors"]) / a["errors"]:+.1f}%), '
                  f'recall {a["recall"]:.4f} -> {b["recall"]:.4f}, precision {a["precision"]:.4f} -> {b["precision"]:.4f}')

print('== HG001: cohort vs alone, same setting (errors PASS, BB)')
for p in sorted({p for g, p, c, s in seen if g == 'HG001'}):
    for s in ('pipeline', 'BP 99.0/99.0', 'BP 99.9/99.0', 'BP 99.9/99.9'):
        for t in ('SNP', 'INDEL'):
            a, b = get('HG001', p, 'single', s, 'BB', t), get('HG001', p, 'cohort', s, 'BB', t)
            if a and b:
                print(f'HG001 {p} [{s}] {t}: alone {a["errors"]:,} (F1 {a["F1"]:.4f}), cohort {b["errors"]:,} (F1 {b["F1"]:.4f}), '
                      f'{100 * (b["errors"] - a["errors"]) / a["errors"]:+.1f}%')

print('== ExcessHet filter: with vs without, HG001 cohort, BP 99.9/99.9, PASS')
for p in sorted({p for g, p, c, s in seen if s == 'BP 99.9/99.9 no ExcessHet'}):
    for region in ('BB', 'LM', 'AD', 'OD'):
        for t in ('SNP', 'INDEL'):
            a, b = get('HG001', p, 'cohort', 'BP 99.9/99.9 no ExcessHet', region, t), get('HG001', p, 'cohort', 'BP 99.9/99.9', region, t)
            if a and b:
                print(f'HG001 {p} {region} {t}: without {a["errors"]:,} (FP {a["FP"]:,}, FN {a["FN"]:,}), with {b["errors"]:,} '
                      f'(FP {b["FP"]:,}, FN {b["FN"]:,}): {b["errors"] - a["errors"]:+,} errors with the filter')

print('== the pipeline VQSR (settings of the joint-calling pipeline) vs Best Practices at the chosen levels 99.9/99.9 (errors PASS, BB and LM)')
for g, p, c in sorted({(g, p, c) for g, p, c, s in seen if s == 'pipeline'}):
    for region in ('BB', 'LM'):
        for t in ('SNP', 'INDEL'):
            a, b = get(g, p, c, 'pipeline', region, t), get(g, p, 'cohort', 'BP 99.9/99.9', region, t)
            if a and b:
                print(f'{g} {p} {c} {region} {t}: pipeline {a["errors"]:,} (F1 {a["F1"]:.4f}) vs BP 99.9/99.9 (cohort) {b["errors"]:,} '
                      f'(F1 {b["F1"]:.4f})')
print('== thresholds by region (calls joint with the cohort, PASS): F1 at a lower threshold minus F1 at 99.9, in points; the SNP threshold at INDEL '
      '99.0 for SNVs and the INDEL threshold at SNP 99.9 for indels (the panel a of Figure 7)')
LEVELS = {'SNP': [('99.0', 'BP 99.0/99.0'), ('99.7', 'BP 99.7/99.0'), ('99.9', 'BP 99.9/99.0')],
          'INDEL': [('99.0', 'BP 99.9/99.0'), ('99.7', 'BP 99.9/99.7'), ('99.9', 'BP 99.9/99.9')]}
REGIONS = ['BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']
for t in ('SNP', 'INDEL'):
    for g, p in (('HG002', 'Illumina'), ('HG002', 'Aviti'), ('HG001', 'Illumina'), ('HG001', 'Aviti')):
        ref = LEVELS[t][2][1]
        for lev, set_ in LEVELS[t][:2]:
            d = {}
            for region in REGIONS:
                a, b = get(g, p, 'cohort', set_, region, t), get(g, p, 'cohort', ref, region, t)
                if a and b:
                    d[region] = 100 * (a['F1'] - b['F1'])
            if not d:
                print(f'{g} {p} {t} {lev}: not run')
                continue
            better = ', '.join(f'{r} {v:+.1f}' for r, v in d.items() if v > 0.5)
            print(f'{g} {p} {t} {lev} minus 99.9: BB {d["BB"]:+.1f}, LM {d["LM"]:+.1f}, AD {d["AD"]:+.1f}; lowest {min(d.values()):+.1f} '
                  f'({min(d, key=d.get)}); higher than 99.9 by more than 0.5 points: {better or "none"}')
print('== SNV rows of Figure 7 use INDEL 99.0 (pair 99.9/99.0 drawn as "99.9"): SNV results there against INDEL 99.9 (the chosen pair 99.9/99.9), PASS')
worst = 0.0
for g, p, c in sorted({(g, p, c) for g, p, c, s in seen if s == 'BP 99.9/99.9'}):
    for region in ('BB', 'LM', 'AD'):
        a, b = get(g, p, c, 'BP 99.9/99.0', region, 'SNP'), get(g, p, c, 'BP 99.9/99.9', region, 'SNP')
        if a and b:
            d = 100 * (b['errors'] - a['errors']) / a['errors']
            worst = max(worst, abs(d))
            print(f'{g} {p} {c} {region} SNP: errors {a["errors"]:,} (INDEL 99.0) -> {b["errors"]:,} (INDEL 99.9), {d:+.2f}%; F1 {a["F1"]:.4f} -> {b["F1"]:.4f}')
print(f'largest change of the SNV errors between INDEL 99.0 and INDEL 99.9: {worst:.2f}%')
print('== the level pair with the highest F1 (PASS) in BB and LM, per series: the five pairs run on every genome and platform')
PAIRS = ['BP 99.0/99.0', 'BP 99.7/99.0', 'BP 99.9/99.0', 'BP 99.9/99.7', 'BP 99.9/99.9']
for g, p, c in sorted({(g, p, c) for g, p, c, s in seen if s == 'BP 99.9/99.9'}):
    for t in ('SNP', 'INDEL', 'ALL'):
        for region in ('BB', 'LM'):
            f = {s: get(g, p, c, s, region, t)['F1'] for s in PAIRS if get(g, p, c, s, region, t)}
            best = max(f, key=f.get)
            chosen = f['BP 99.9/99.9']
            others = ', '.join(f'{s[3:]} {v:.4f}' for s, v in f.items() if s != 'BP 99.9/99.9')
            print(f'{g} {p} {c} {region} {t}: best {best[3:]} F1 {f[best]:.4f}; 99.9/99.9 {chosen:.4f} (difference {100 * (chosen - f[best]):+.2f} points); others {others}')
print('F7_TABLE_DONE')
