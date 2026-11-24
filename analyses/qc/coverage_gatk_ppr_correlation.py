#!/usr/bin/env python3
"""Pearson correlation of the mean coverage of the GATK and the PPR alignments over the individuals of each platform
(Supplementary Note 1). Excluded samples are left out; Aviti uses the GATK table numbered by individual.
Usage: python3 coverage_gatk_ppr_correlation.py [qc folder with tables_GATKbam_qc/ and tables_PPRbam_qc/]
"""
import csv
import os
import re
import sys
from statistics import correlation

QC = sys.argv[1] if len(sys.argv) > 1 else '.'


def individual(name):
    """the number of the individual in a QC sample name (its first digits; numbered names pass unchanged)"""
    return re.search(r'\d+', name).group()


def mean_coverage(path):
    with open(path, newline='') as fh:
        return {individual(r['sample']): float(r['mean_coverage']) for r in csv.DictReader(fh, delimiter='\t') if r['excluded'] == 'no'}


for plat in ('illumina', 'aviti'):
    g = mean_coverage(os.path.join(QC, 'tables_GATKbam_qc', f'bam_qc_per_sample_{plat}_gatk.tsv'))
    p = mean_coverage(os.path.join(QC, 'tables_PPRbam_qc', f'bam_qc_per_sample_{plat}_ppr.tsv'))
    both = sorted(set(g) & set(p))
    print(f'{plat}: n = {len(both)} (GATK {len(g)}, PPR {len(p)}), Pearson r = {correlation([g[s] for s in both], [p[s] for s in both]):.4f}')
