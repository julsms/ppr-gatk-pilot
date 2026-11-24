#!/usr/bin/env python3
"""The graph decomposition on HG002 in one set of tables: total errors (SNP FN + SNP FP + INDEL FN + INDEL FP, hap.py PASS rows, as
Siren et al. 2024) and F1 of every set-up, the effect of each step (aligner, graph, filtering, personalization) and the share of each
step in the change from BWA-MEM to PPR. PLATFORM=illumina (default) or aviti.
Usage: python3 decompose.py   (in the benchmark folder, or BENCHMARK_REVISION_DIR=<folder>; README.md has the layout)
"""
import csv
import os
import sys
from collections import defaultdict

BR = os.environ.get('BENCHMARK_REVISION_DIR', '.')
PLATFORM = os.environ.get('PLATFORM') or 'illumina'      # unset or empty = illumina
if PLATFORM not in ('illumina', 'aviti'):
    sys.exit(f"PLATFORM must be illumina or aviti, not '{PLATFORM}'")
PFX = {'illumina': 'Illumina', 'aviti': 'Aviti'}[PLATFORM]   # prefix of the arm folder names
D = os.path.join(BR, 'comparisons', f'Decomposition_summary_with_Siren_{PFX}')   # ..._Illumina or ..._Aviti
ARMS_TSV = os.path.join(BR, 'scripts', 'graph_arms', 'arms.tsv')
A0 = f'{PFX}_GATK_Deepvariant_markdup'
A1 = f'{PFX}_giraffe_GRCh38only_Deepvariant'
A2 = f'{PFX}_giraffe_HPRC_d9_Deepvariant'
A2B = f'{PFX}_giraffe_HPRC_default_Deepvariant'
A3 = f'{PFX}_PPR_Deepvariant'
A3GBZ = f'{PFX}_PPR_surjectGBZ_Deepvariant'
GATK = f'{PFX}_GATK_HaplotypeCaller'
# short labels in the output tables
CODE = {A0: 'BWA-MEM', A1: 'giraffe_GRCh38only', A2: 'giraffe_HPRC_d9', A2B: 'giraffe_HPRC_default', A3: 'PPR', A3GBZ: 'PPR_surjectGBZ',
        GATK: 'GATK_HC'}
# other names of the Illumina arms, read as their folder names (arms.tsv uses the folder names)
OLD = {'bwa_dv_matched': A0, 'bwa_dv_default': A0, 'giraffe_grch38': A1, 'giraffe_hprc_d9': A2, 'giraffe_hprc_default': A2B,
       'ppr_sampled': A3, 'ppr_sampled_gbz_route': A3GBZ, 'gatk_hc': GATK} if PLATFORM == 'illumina' else {}
METHODS = ['ga4gh_v5.0q', 'ga4gh_cmrg_v1.00']
WHOLE = {'ga4gh_v5.0q': 'BB', 'ga4gh_cmrg_v1.00': 'CMRG_v1.00_truth'}
TYPES = ['ALL', 'SNP', 'INDEL']
EFFECTS = [   # name, from, to, what changes
    ('caller', GATK, A0, 'HaplotypeCaller -> DeepVariant, both on the BWA-MEM BAM'),
    ('aligner', A0, A1, 'BWA-MEM -> giraffe, both on linear GRCh38 (and duplicate-flagged reads ignored -> no duplicate flags)'),
    ('graph_d9', A1, A2, 'GRCh38-only graph -> HPRC v1.1 d9 graph, no sampling'),
    ('graph_default', A1, A2B, 'GRCh38-only graph -> HPRC v1.1 default graph, no sampling'),
    ('filtering', A2B, A2, 'default graph -> d9 graph'),
    ('personalization_vs_d9', A2, A3, 'd9 graph -> diploid sampled graph'),
    ('personalization_vs_default', A2B, A3, 'default (full) graph -> diploid sampled graph'),
    ('bwa_to_d9', A0, A2, 'BWA-MEM + DeepVariant -> giraffe d9 (aligner and graph together)'),
    ('total_vs_bwa', A0, A3, 'BWA-MEM + DeepVariant -> PPR'),
    ('total_vs_gatk', GATK, A3, 'GATK best practice -> PPR'),
]
CHAINS = {'via_d9': [A0, A1, A2, A3], 'via_default': [A0, A1, A2B, A3]}
NAN = float('nan')


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def fl(v):
    return float(v) if v not in (None, '') else NAN


def ni(v):
    return int(float(v)) if v not in (None, '') else 0


def f1(tp_t, fn, tp_q, fp):
    rec = tp_t / (tp_t + fn) if tp_t + fn else NAN
    pre = tp_q / (tp_q + fp) if tp_q + fp else NAN
    return rec, pre, (2 * rec * pre / (rec + pre) if rec + pre else NAN)


def pct(a, b):
    return 100.0 * (b - a) / a if a else NAN


def fmt(x, d=2):
    return 'NA' if x is None or x != x else f'{x:.{d}f}'


def sfmt(x, d=6):
    return 'NA' if x is None or x != x else f'{x:+.{d}f}'


def label(a):
    return OLD.get(a, a)


def code(a):
    return CODE.get(a, a)


arms, desc = [], {}
for r in read_tsv(ARMS_TSV):
    a = label(r['arm'])
    if not a.startswith(PFX + '_'):     # the rows of the platform
        continue
    if a not in arms:
        arms.append(a)
        desc[a] = r['description']
os.makedirs(D, exist_ok=True)

route_note, alias = '', {}
rc = os.path.join(BR, 'checks', 'surject_GBZ_vs_XG', 'results', 'route_check.txt')
if PLATFORM == 'illumina' and os.path.exists(rc):   # the route check was made on Illumina only
    route = open(rc).read().strip()
    route_note = f'Surjection route check: {route}'
    if route.startswith('ROUTE_DIFFERENT'):
        if os.path.exists(os.path.join(BR, A3GBZ, 'BENCH_DONE')):
            alias[A3] = A3GBZ
            route_note += f' -> the PPR call set surjected from the GBZ ({A3GBZ}) is used in place of {A3}.'
        else:
            route_note += (f' -> WARNING: {A3GBZ} has no finished benchmark yet; {A3} here was surjected by another route than the graph '
                           'arms.')


def tables(arm):
    """(method, row) of the finished benchmark tables of one call set"""
    out = []
    g = os.path.join(BR, arm, 'summary.tsv')
    if os.path.exists(os.path.join(BR, arm, 'BENCH_DONE')) and os.path.exists(g):
        for x in read_tsv(g):
            out.append(('ga4gh_cmrg_v1.00' if x['region'] == 'CMRG_v1.00_truth' else 'ga4gh_v5.0q', x))
    return out


# metrics[(arm, method, region, type)] = dict
metrics, warnings = {}, []
for arm in arms:
    acc = defaultdict(lambda: [0, 0, 0, 0])
    for m, x in tables(arm):
        if x['filter'] != 'PASS' or x['type'] not in ('SNP', 'INDEL'):
            continue
        tp_t, fn, fp = ni(x['TP']), ni(x['FN']), ni(x['FP'])
        tp_q = ni(x['query_total']) - fp - ni(x['query_unk'])
        prec = fl(x['precision'])
        if tp_q + fp and prec == prec and abs(tp_q / (tp_q + fp) - prec) > 1e-4:
            warnings.append(f"{arm} {m} {x['region']} {x['type']}: precision from the counts {tp_q / (tp_q + fp):.6f}, hap.py {prec:.6f}")
        metrics[(arm, m, x['region'], x['type'])] = {'truth_TP': tp_t, 'FN': fn, 'query_TP': tp_q, 'FP': fp, 'errors': fn + fp,
                                                    'recall': fl(x['recall']), 'precision': prec, 'F1': fl(x['F1'])}
        s = acc[(arm, m, x['region'])]
        s[0] += tp_t; s[1] += fn; s[2] += tp_q; s[3] += fp
    for (a, m, reg), (tp_t, fn, tp_q, fp) in acc.items():
        rec, pre, f = f1(tp_t, fn, tp_q, fp)
        metrics[(a, m, reg, 'ALL')] = {'truth_TP': tp_t, 'FN': fn, 'query_TP': tp_q, 'FP': fp, 'errors': fn + fp, 'recall': rec,
                                       'precision': pre, 'F1': f}
present = sorted({k[0] for k in metrics})


def get(arm, m, reg, t):
    return metrics.get((alias.get(arm, arm), m, reg, t))


regions = defaultdict(list)
for (a, m, reg, t) in metrics:
    if reg not in regions[m]:
        regions[m].append(reg)

with open(os.path.join(D, 'arm_metrics.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['label', 'arm', 'method', 'region', 'type', 'truth_TP', 'FN', 'query_TP', 'FP', 'errors', 'recall', 'precision', 'F1'])
    for k in sorted(metrics):
        v = metrics[k]
        w.writerow([code(k[0])] + list(k) + [v['truth_TP'], v['FN'], v['query_TP'], v['FP'], v['errors'], fmt(v['recall'], 6),
                                             fmt(v['precision'], 6), fmt(v['F1'], 6)])

EFF_HEAD = ['effect', 'from_label', 'to_label', 'from_arm', 'to_arm', 'method', 'region', 'type', 'errors_from', 'errors_to', 'delta_errors',
            'rel_change_pct', 'delta_FN', 'delta_FP', 'F1_from', 'F1_to', 'delta_F1']
eff_rows = []
for name, a, b, _ in EFFECTS:
    for m in METHODS:
        for reg in regions.get(m, []):
            for t in TYPES:
                x, y = get(a, m, reg, t), get(b, m, reg, t)
                if not x or not y:
                    continue
                eff_rows.append(dict(zip(EFF_HEAD, [name, code(a), code(b), alias.get(a, a), alias.get(b, b), m, reg, t, x['errors'],
                                                    y['errors'], y['errors'] - x['errors'], fmt(pct(x['errors'], y['errors'])),
                                                    y['FN'] - x['FN'], y['FP'] - x['FP'], fmt(x['F1'], 6), fmt(y['F1'], 6),
                                                    sfmt(y['F1'] - x['F1'], 6)])))
with open(os.path.join(D, 'effects.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=EFF_HEAD, delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(eff_rows)

SHARE_HEAD = ['chain', 'method', 'region', 'type', 'errors_BWA-MEM', 'errors_PPR', 'total_delta', 'total_rel_change_pct', 'aligner_delta',
              'aligner_share_pct', 'graph_delta', 'graph_share_pct', 'personalization_delta', 'personalization_share_pct']
share_rows = []
for chain, arms_c in CHAINS.items():
    for m in METHODS:
        for reg in regions.get(m, []):
            for t in TYPES:
                vals = [get(a, m, reg, t) for a in arms_c]
                if not all(vals):
                    continue
                e = [v['errors'] for v in vals]
                total = e[-1] - e[0]
                steps = [e[i + 1] - e[i] for i in range(3)]
                assert sum(steps) == total
                share_rows.append([chain, m, reg, t, e[0], e[-1], total, fmt(pct(e[0], e[-1]))] +
                                  [x for s in steps for x in (s, fmt(100.0 * s / total if total else NAN, 1))])
with open(os.path.join(D, 'shares.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(SHARE_HEAD)
    w.writerows(share_rows)

SIDE_HEAD = ['effect', 'from -> to', 'what changes', 'v5.0q GA4GH (%)', 'CMRG v1.00 (%)']
side = []
for name, a, b, text in EFFECTS:
    row = [name, f'{code(a)} -> {code(b)}', text]
    for m in METHODS:
        x, y = get(a, m, WHOLE[m], 'ALL'), get(b, m, WHOLE[m], 'ALL')
        row.append(fmt(pct(x['errors'], y['errors'])) if x and y else 'NA')
    side.append(row)

# readable report
TITLE = {'illumina': 'Illumina NovaSeq PCR-free SRR14724532, trimmed', 'aviti': 'Element Aviti, Adept library SRR33551004, trimmed'}
L = [f'# Decomposition on HG002 ({TITLE[PLATFORM]})', '',
     'Total errors = SNP FN + SNP FP + INDEL FN + INDEL FP (hap.py PASS rows), as in Siren et al. 2024 Supplementary Table 5. '
     'Relative change = (errors of the second arm - errors of the first) / errors of the first. Benchmarks: GA4GH way on GIAB v5.0q '
     '(happy.sh: BB and the strata of strata_revision; CMRG v1.00 with its own truth).', '']
if route_note:
    L += [route_note, '']
L += ['## Arms (whole benchmark region, PASS, SNP + INDEL)', '',
      '| label | arm | errors v5.0q (GA4GH) | errors CMRG v1.00 | F1 SNP v5.0q | F1 INDEL v5.0q | description |',
      '|---|---|---|---|---|---|---|']
for a in arms:
    if a not in present:
        continue
    y, c = metrics.get((a, 'ga4gh_v5.0q', 'BB', 'ALL')), metrics.get((a, 'ga4gh_cmrg_v1.00', 'CMRG_v1.00_truth', 'ALL'))
    s, i = metrics.get((a, 'ga4gh_v5.0q', 'BB', 'SNP')), metrics.get((a, 'ga4gh_v5.0q', 'BB', 'INDEL'))
    L.append(f"| {code(a)} | {a} | {y['errors'] if y else 'NA'} | {c['errors'] if c else 'NA'} | "
             f"{fmt(s['F1'], 6) if s else 'NA'} | {fmt(i['F1'], 6) if i else 'NA'} | {desc.get(a, '')} |")
missing = [a for a in arms if a not in present]
if missing:
    L += ['', 'Not benchmarked (yet): ' + ', '.join(missing)]
L += ['', '## Effects (whole benchmark region, SNP + INDEL, relative change in total errors, %)', '',
      '| ' + ' | '.join(SIDE_HEAD) + ' |', '|---' * len(SIDE_HEAD) + '|']
for row in side:
    L.append('| ' + ' | '.join(row) + ' |')
L += ['', '## Shares of the change from BWA-MEM + DeepVariant to PPR (whole region, SNP + INDEL)', '',
      '| chain | method | errors BWA-MEM | errors PPR | total change | aligner | graph | personalization |', '|---|---|---|---|---|---|---|---|']
for r in share_rows:
    if r[2] in WHOLE.values() and r[3] == 'ALL':
        L.append(f'| {r[0]} | {r[1]} | {r[4]} | {r[5]} | {r[6]} ({r[7]} %) | {r[8]} ({r[9]} %) | {r[10]} ({r[11]} %) | {r[12]} ({r[13]} %) |')
L += ['', '## Main effects per region (GA4GH way, PASS; relative change in total errors, %; F1 change in brackets)', '']
main = ['aligner', 'graph_d9', 'graph_default', 'personalization_vs_d9', 'personalization_vs_default']
idx = {(r['effect'], r['method'], r['region'], r['type']): r for r in eff_rows}
reg_list = [('ga4gh_v5.0q', reg) for reg in regions.get('ga4gh_v5.0q', [])] + \
           [('ga4gh_cmrg_v1.00', reg) for reg in regions.get('ga4gh_cmrg_v1.00', [])]
for t in TYPES:
    L += [f'### {t}', '', '| region | ' + ' | '.join(main) + ' |', '|---' * (len(main) + 1) + '|']
    for m, reg in reg_list:
        cells = []
        for e in main:
            r = idx.get((e, m, reg, t))
            cells.append(f"{r['rel_change_pct']} ({r['delta_F1']})" if r else 'NA')
        L.append(f'| {reg} | ' + ' | '.join(cells) + ' |')
    L.append('')
if warnings:
    L += ['## Check of the counts', '', 'Query TP computed as query_total - FP - query_unk gives a precision that differs from hap.py\'s '
          'by more than 0.0001 in these rows (the F1 of type ALL depends on it):', ''] + [f'- {x}' for x in warnings] + ['']
L += ['## Caveats', '',
      '- BWA-MEM -> giraffe_GRCh38only (aligner) also contains duplicate handling: the BWA-MEM BAM is duplicate-marked and DeepVariant '
      'with the PPR settings (keep_duplicates not set) ignores duplicate-flagged reads, while the giraffe BAMs carry no duplicate flags. '
      'The BWA-MEM call set is scored with the Ensembl-named truth, bed and strata of happy.sh, and renamed to chr names for CMRG '
      'v1.00.',
      {'illumina': '- HG002 is in the training data of the DeepVariant WGS model (1.9); one sample; our reads are about 35x before '
                   'trimming.',
       'aviti': '- HG002 is in the training data of the DeepVariant WGS model (1.9); one sample; our Aviti reads are about 46x '
                'before trimming (142.1 Gb).'}[PLATFORM],
      '- The graphs contain GRCh38 from the no-alt analysis set (chrY PAR masked) while BWA-MEM used the Ensembl primary assembly; PPR '
      'does not run the vg_wdl post-processing (left-alignment, indel realignment).']
if PLATFORM == 'aviti':
    L += ['- The PPR call set mapped with 24 threads; the giraffe arms mapped with 64 threads on other machines. Thread count and '
          'host are run settings, not a change of graph or aligner, but the giraffe run times are not comparable between any of these '
          'set-ups (comparisons/Mapping_metrics_graph_arms_Aviti).',
          '- No GATK HaplotypeCaller arm for Aviti in this decomposition, so the caller and total_vs_gatk rows are NA by design (not a '
          'failed benchmark). The Aviti GATK HaplotypeCaller call set is compared with BWA-MEM + DeepVariant and PPR in the pipeline '
          'benchmark.',
          '- The surjection route check found the same records for the two routes on the Illumina PPR alignments; it was '
          'not repeated for Aviti, whose PPR VCF is used as delivered.']
with open(os.path.join(D, 'decomposition.md'), 'w') as fh:
    fh.write('\n'.join(L) + '\n')
print('\n'.join(L))
