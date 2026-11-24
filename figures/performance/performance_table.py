#!/usr/bin/env python3
"""Data of Supplementary Figure S28 and Supplementary Note 2: execution time, CPU time and memory per sample of GATK and PPR (samples
numbered), and their n, mean, SD and range. GATK runs with an incomplete trace (gatk_run_completeness.tsv) are left out of every metric.
Inputs: <TIMELINE_DIR> (analyses/performance): the three .csv tables, gatk_run_completeness.tsv, ppr_v3_matched_gatk_summary.tsv.
Outputs: performance_data.tsv, performance_summary.tsv and performance_data.provenance.txt next to this script.
"""
import csv
import hashlib
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TD = os.environ.get('TIMELINE_DIR', '.')
# metric: (file, column, unit factor, decimals of the AVERAGE row)
METRICS = {'time_h': ('total_time_results.csv', 'Total_Time_Hours', 1.0, 2),
           'cpu_h': ('cpu_minutes_results.csv', 'CPU_Minutes', 1 / 60, 1),
           'max_rss_gb': ('memory_usage_results.csv', 'Max_Peak_RSS_GB', 1.0, 1),
           'total_rss_gb': ('memory_usage_results.csv', 'Total_Peak_RSS_GB', 1.0, 1)}
PIPELINES = [('GATK', 'GATK'), ('PPR', 'PPR')]
PLATFORMS = [('illumina', 'Illumina'), ('aviti', 'Aviti')]


def md5(path):
    with open(path, 'rb') as fh:
        return hashlib.md5(fh.read()).hexdigest()


def fail(msg):
    print('ERROR:', msg)
    print('PERFORMANCE_TABLE_FAILED')
    sys.exit(1)


data, summary, inputs = [], [], []
path = os.path.join(TD, 'gatk_run_completeness.tsv')
inputs.append(path)
complete = {'illumina': set(), 'aviti': set()}
incomplete = {'illumina': set(), 'aviti': set()}
with open(path, newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        (complete if r['complete'] == 'yes' else incomplete)[r['platform']].add(r['sample'])
for metric, (fname, col, factor, dec) in METRICS.items():
    path = os.path.join(TD, fname)
    if path not in inputs:
        inputs.append(path)
    with open(path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    for pipe, plabel in PIPELINES:
        for plat, tlabel in PLATFORMS:
            key = f'{plabel}_{tlabel}'
            vals = {r['Sample']: float(r[col]) for r in rows if r['Pipeline'] == key and r['Sample'] != 'AVERAGE'}
            avg = [float(r[col]) for r in rows if r['Pipeline'] == key and r['Sample'] == 'AVERAGE']
            if not vals or len(avg) != 1:
                fail(f'{fname} {key}: no samples or no single AVERAGE row')
            mean_raw = statistics.mean(vals.values())
            if abs(round(mean_raw, dec) - avg[0]) > 10 ** -dec * 1.01:
                fail(f'{fname} {key}: mean {mean_raw:.4f} is not the AVERAGE row {avg[0]}')
            if pipe == 'GATK':
                vals = {s: x for s, x in vals.items() if s not in incomplete[plat]}
                if set(vals) != complete[plat]:
                    fail(f'{fname} {key}: the samples are not the complete GATK runs of gatk_run_completeness.tsv')
            v = [x * factor for _, x in sorted(vals.items())]
            for i, x in enumerate(v, 1):
                data.append([pipe, plat, i, metric, f'{x:.4f}'])
            summary.append([pipe, plat, metric, len(v), f'{statistics.mean(v):.4f}', f'{statistics.stdev(v):.4f}', f'{min(v):.4f}',
                            f'{max(v):.4f}'])
path = os.path.join(TD, 'ppr_v3_matched_gatk_summary.tsv')
inputs.append(path)
with open(path, newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        # only the execution-time pairs: GATK has a time for its complete runs only, while the CPU and memory pairs of this table
        # include the incomplete GATK runs
        if r['metric'] == 'time':
            summary.append(['same samples', r['platform'], r['metric'], r['n_pairs'], f'PPR {r["mean_PPR_v3"]}',
                            f'GATK {r["mean_GATK"]}', f'ratio {r["ratio_PPR_over_GATK"]}', ''])
outs = []
for name, head, rows in (('performance_data.tsv', ['pipeline', 'platform', 'sample_index', 'metric', 'value'], data),
                         ('performance_summary.tsv', ['pipeline', 'platform', 'metric', 'n', 'mean', 'sd', 'min', 'max'], summary)):
    out = os.path.join(HERE, name)
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(head)
        w.writerows(rows)
    outs.append(out)
    print(f'{len(rows)} rows written to {out}')
with open(os.path.join(HERE, 'performance_data.provenance.txt'), 'w') as fh:
    for p in inputs:
        fh.write(f'input\t{md5(p)}\t{p}\n')
    fh.write(f'script\t{md5(os.path.abspath(__file__))}\t{os.path.abspath(__file__)}\n')
    for plat in ('illumina', 'aviti'):
        fh.write(f'excluded\tGATK runs with an incomplete trace, {plat}\t{len(incomplete[plat])}\n')
    for p in outs:
        fh.write(f'table\t{md5(p)}\t{p}\n')
print('PERFORMANCE_TABLE_OK')
