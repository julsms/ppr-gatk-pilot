#!/usr/bin/env python3
"""Alignment metrics of each set-up: from vg stats -a on the GAM (aligned, perfect, gapless, properly paired, MAPQ 60) and from
bam_metrics.sh on the BAM (mapped, properly paired, MAPQ >= 60, exact, gapless), with the giraffe time and memory from the giraffe log.
The paths come from scripts/graph_arms/arms.tsv. PLATFORM=illumina (default) or aviti.
Usage: python3 mapping_summary.py   (in the benchmark folder, or BENCHMARK_REVISION_DIR=<folder>)
"""
import csv
import os
import re
import sys

RUN_BR = '.'
BR = os.environ.get('BENCHMARK_REVISION_DIR', RUN_BR)
PLATFORM = os.environ.get('PLATFORM') or 'illumina'      # unset or empty = illumina
if PLATFORM not in ('illumina', 'aviti'):
    sys.exit(f"PLATFORM must be illumina or aviti, not '{PLATFORM}'")
PFX = {'illumina': 'Illumina', 'aviti': 'Aviti'}[PLATFORM]   # prefix of the arm folder names
D = os.path.join(BR, 'comparisons', f'Mapping_metrics_graph_arms_{PFX}')   # ..._Illumina or ..._Aviti
ARMS_TSV = os.path.join(BR, 'scripts', 'graph_arms', 'arms.tsv')
TITLE = {'illumina': 'HG002 Illumina PCR-free trimmed', 'aviti': 'HG002 Aviti SRR33551004 trimmed'}


def at_br(p):
    """arms.tsv paths start with ./ (the benchmark folder); with BENCHMARK_REVISION_DIR set they are read from that folder"""
    if p and BR != RUN_BR and p.startswith(RUN_BR + '/'):
        return BR + p[len(RUN_BR):]
    return p


def vg_stats(p):
    p = at_br(p)
    if not p or p == '-' or not os.path.exists(p):
        return {}
    out = {}
    for line in open(p):
        k, _, v = line.partition(':')
        v = v.strip()
        if k in ('Total primary', 'Total aligned', 'Total perfect', 'Total gapless (softclips allowed)', 'Total properly paired'):
            out[k] = int(v.split()[0])
        elif k == 'Mapping quality':
            m = re.search(r'max (\d+) \((\d+) reads\)', v)
            if m and m.group(1) == '60':
                out['mapq60'] = int(m.group(2))
    return out


def align_err(p):
    p = at_br(p)
    if not p or p == '-' or not os.path.exists(p):
        return {}
    out = {}
    for line in open(p, errors='replace'):
        m = re.search(r'Mapped (\d+) reads across (\d+) threads in ([0-9.]+) seconds', line)
        if m:
            out['map_seconds'] = float(m.group(3))
            out['map_threads'] = int(m.group(2))
        m = re.search(r'Memory footprint: ([0-9.]+) GB', line)
        if m:
            out['memory_gb'] = float(m.group(1))
        m = re.search(r'Loading and initialization: ([0-9.]+) seconds', line)
        if m:
            out['load_seconds'] = float(m.group(1))
    return out


def bam(arm):
    p = os.path.join(D, 'bam', f'{arm}.tsv')
    if not os.path.exists(p):
        return {}
    return {k: int(v) for k, v in (line.rstrip('\n').split('\t') for line in open(p) if line.strip())}


def flagstat(p):
    """primary records, primary mapped and properly paired from samtools flagstat (1.13 or later), QC-passed + QC-failed, as
    bam_metrics.sh counts them; {} without a primary line"""
    p = at_br(p)
    if not p or p == '-' or not os.path.exists(p):
        return {}
    out = {}
    for line in open(p):
        m = re.match(r'(\d+) \+ (\d+) (primary|primary mapped|properly paired)( \(|$)', line.strip())
        if m:
            out[{'primary': 'primary', 'primary mapped': 'mapped', 'properly paired': 'properly_paired'}[m.group(3)]] = \
                int(m.group(1)) + int(m.group(2))
    return out if 'primary' in out else {}


def pc(a, b):
    return f'{100.0 * a / b:.2f}' if a is not None and b else 'NA'


HEAD = ['arm', 'graph_primary', 'graph_aligned_pct', 'graph_perfect_pct', 'graph_gapless_pct', 'graph_properly_paired_pct',
        'graph_mapq60_pct', 'bam_primary', 'bam_mapped_pct', 'bam_properly_paired_pct', 'bam_mapq60_pct', 'bam_exact_pct',
        'bam_gapless_pct', 'giraffe_map_seconds', 'giraffe_threads', 'giraffe_load_seconds', 'giraffe_memory_gb']
if PLATFORM == 'aviti':
    HEAD = HEAD + ['giraffe_host']     # Aviti only


def staged_check(p):
    """the check column of the last row of Aviti_PPR_mapping_stats/run_info/source.tsv for the file p, or ''"""
    src = os.path.join(BR, 'Aviti_PPR_mapping_stats', 'run_info', 'source.tsv')
    if not p or p == '-' or not os.path.exists(src):
        return ''
    out = ''
    with open(src) as fh:
        for x in csv.DictReader(fh, delimiter='\t'):
            if x.get('file') == at_br(p):
                out = x.get('check') or ''
    return out


rows = [HEAD]
from_flagstat, markdup_flagstat = [], []
with open(ARMS_TSV) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if not r['arm'].startswith(PFX + '_'):   # the rows of the platform
            continue
        g, b, e = vg_stats(r['vg_stats']), bam(r['arm']), align_err(r['align_err'])
        if not b:
            b = flagstat(r.get('flagstat'))      # no BAM metrics: mapped and properly paired from flagstat, the rest NA
            if b:
                from_flagstat.append(r['arm'])
                if PLATFORM == 'aviti' and staged_check(r.get('flagstat')) == 'markdup_flagstat_fallback':
                    markdup_flagstat.append(r['arm'])
        if not (g or b or e):
            continue
        gp, bp = g.get('Total primary'), b.get('primary')
        if not b:
            exact = 'NA'
        elif b.get('no_NM_tag', 0) == 0:
            exact = pc(b.get('exact'), bp)
        else:
            exact = f"NA ({b['no_NM_tag']} records without NM tag)"
        rows.append([r['arm'], gp if gp else 'NA', pc(g.get('Total aligned'), gp), pc(g.get('Total perfect'), gp),
                     pc(g.get('Total gapless (softclips allowed)'), gp), pc(g.get('Total properly paired'), gp), pc(g.get('mapq60'), gp),
                     bp if bp else 'NA', pc(b.get('mapped'), bp), pc(b.get('properly_paired'), bp), pc(b.get('mapq60'), bp), exact,
                     pc(b.get('gapless'), bp), e.get('map_seconds', 'NA'), e.get('map_threads', 'NA'), e.get('load_seconds', 'NA'),
                     e.get('memory_gb', 'NA')] + ([r.get('giraffe_host') or 'NA'] if PLATFORM == 'aviti' else []))
os.makedirs(D, exist_ok=True)
with open(os.path.join(D, 'mapping_metrics.tsv'), 'w', newline='') as fh:
    csv.writer(fh, delimiter='\t', lineterminator='\n').writerows(rows)

L = [f'# Mapping metrics per arm ({TITLE[PLATFORM]}; context for the decomposition)', '',
     'Graph space = vg stats -a on the GAM (percent of primary reads); BAM = the surjected or BWA-MEM BAM (percent of primary records).', '']
notes = []
if from_flagstat:
    notes.append('BAM columns of ' + ', '.join(from_flagstat) + ': from samtools flagstat of the BAM (no BAM metrics made): primary '
                 'records, mapped and properly paired only; MAPQ >= 60, exact and gapless are NA.')
if markdup_flagstat:
    notes.append('The flagstat of ' + ', '.join(markdup_flagstat) + ' is the one of the same BAM after Picard MarkDuplicates: '
                 'duplicates are flagged, not removed, so primary records, mapped and properly paired are those of the BAM.')
threads = sorted({str(x[HEAD.index('giraffe_threads')]) for x in rows[1:]} - {'NA'})
hosts = sorted({str(x[HEAD.index('giraffe_host')]) for x in rows[1:]} - {'NA', '-', ''}) if 'giraffe_host' in HEAD else []
if len(hosts) > 1:
    notes.append('giraffe ran ' + ('with different thread counts (' + ', '.join(threads) + ') and ' if len(threads) > 1 else '') +
                 'on different hosts (' + ', '.join(hosts) + '; columns giraffe_threads, giraffe_host): its mapping times are not '
                 'comparable between these rows.')
elif len(threads) > 1:
    notes.append('giraffe ran with different thread counts (' + ', '.join(threads) + '; column giraffe_threads): its mapping times '
                 'are not comparable between these rows.')
L += [f'- {x}' for x in notes] + ([''] if notes else [])
L += ['| ' + ' | '.join(HEAD) + ' |', '|---' * len(HEAD) + '|']
L += ['| ' + ' | '.join(str(x) for x in r) + ' |' for r in rows[1:]]
with open(os.path.join(D, 'mapping_metrics.md'), 'w') as fh:
    fh.write('\n'.join(L) + '\n')
print('\n'.join(L))
