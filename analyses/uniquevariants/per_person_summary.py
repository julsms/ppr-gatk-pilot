#!/usr/bin/env python3
"""per_person_summary.py - the numbers of the Results sentences on the novel records per individual, from the per_person tables of
novel_per_person.py: per platform, the novel records summed over the 23 individuals, the share held by the three individuals with the
most, the range of those three and the range of the other 20. Prints no individual value.
Inputs : one or more per_person tables (label=path)
Outputs: stdout; last line PER_PERSON_SUMMARY_DONE
Usage  : python3 per_person_summary.py GATKvsPPR=per_person.tsv PPRvsGATK=per_person_PPRvsGATK.tsv PPRvsDV=per_person_PPRvsDV.tsv \
             DVvsPPR=per_person_DVvsPPR.tsv
"""
import csv
import sys

for arg in sys.argv[1:]:
    label, path = arg.split('=', 1)
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    for plat in ('illumina', 'aviti'):
        v = sorted((int(r['novel']) for r in rows if r['platform'] == plat), reverse=True)
        top, rest, total = v[:3], v[3:], sum(v)
        print(f'{label} {plat}: individuals {len(v)}; novel records {total:,}; three with the most hold '
              f'{100 * sum(top) / total:.1f}% ({min(top):,} to {max(top):,}); the other {len(rest)}: {min(rest):,} to {max(rest):,}')
print('PER_PERSON_SUMMARY_DONE')
