#!/bin/bash
# For each unique set: the individuals that carry each allele in the full call sets, the singletons, and the multisample VCF
# (bcftools merge -m none). Run after extract_unique.sh for every individual and platform.
# Usage: bash carriers.sh
set -euo pipefail

export CONC_OUT=${CONC_OUT:-concordance}
export UV_OUT=${UV_OUT:-uniquevariants}
export RELATED=${RELATED:-6}        # one of the two first-degree relatives, left out of the singleton counts
HERE=$(cd "$(dirname "$0")" && pwd)

# 02_carrier_table.py counts the full call set of an individual only when a .DONE file sits next to its key table
for k in $CONC_OUT/prep/*/*.keys.tsv.gz; do touch ${k%.keys.tsv.gz}.DONE; done
for P in illumina aviti; do
    for U in PPRvsGATK GATKvsPPR PPRvsDV DVvsPPR PPRvsBOTH; do
        python3 $HERE/02_carrier_table.py --uset ${U}_$P
    done
done
