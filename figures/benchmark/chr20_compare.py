#!/usr/bin/env python3
"""The chromosome 20-22 numbers of the Results: for each platform and type, F1 of GATK, BWA-MEM + DeepVariant and PPR on the whole
genome and on chromosomes 20-22, the order of the three, and the errors (FN + FP) of DV and PPR relative to GATK (percent fewer
errors), PASS calls. Input: chr20_22_benchmark.tsv (chr20_table.py). Output: a tab-separated table on stdout.
Usage: python3 chr20_compare.py
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PIPES = ['GATK', 'BWA-MEM + DeepVariant', 'PPR']
t = {(r['platform'], r['pipeline'], r['scope'], r['type']): r for r in
     csv.DictReader(open(os.path.join(HERE, 'chr20_22_benchmark.tsv')), delimiter='\t') if r['filter'] == 'PASS'}
print('platform\ttype\tscope\tF1 GATK\tF1 DV\tF1 PPR\torder\terrors GATK\tDV fewer errors %\tPPR fewer errors %\ttruth (TP+FN)')
for plat in ('Illumina', 'Aviti'):
    for typ in ('SNP', 'INDEL'):
        for scope in ('genome', 'chr20_22'):
            f1 = {p: float(t[(plat, p, scope, typ)]['f1']) for p in PIPES}
            err = {p: int(t[(plat, p, scope, typ)]['FN']) + int(t[(plat, p, scope, typ)]['FP']) for p in PIPES}
            truth = int(t[(plat, 'PPR', scope, typ)]['TP']) + int(t[(plat, 'PPR', scope, typ)]['FN'])
            order = ' < '.join(sorted(PIPES, key=lambda p: f1[p])).replace('BWA-MEM + DeepVariant', 'DV')
            red = {p: 100 * (err['GATK'] - err[p]) / err['GATK'] for p in PIPES[1:]}
            print(f'{plat}\t{typ}\t{scope}\t{f1["GATK"]:.4f}\t{f1["BWA-MEM + DeepVariant"]:.4f}\t{f1["PPR"]:.4f}\t{order}\t'
                  f'{err["GATK"]}\t{red["BWA-MEM + DeepVariant"]:.1f}\t{red["PPR"]:.1f}\t{truth}')
print('CHR20_COMPARE_DONE')
