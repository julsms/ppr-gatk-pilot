#!/bin/bash
# The GIAB HG002 v5.0q truth for the precision of the unique calls: the record whose REF differs from the reference is removed,
# then the truth is prepared as a call set (one allele per record, analysis space, key table with the strata). Run once.
# Usage: bash hg002_truth.sh <HG002_GRCh38_v5.0q_smvar.vcf.gz>
set -euo pipefail

TRUTH=$1
export REF=${REF:-GRCh38_chr.fa}    # GRCh38 with chr names (+ .fai), as in the concordance step
H=${UV_OUT:-uniquevariants}/hg002
CONC=$(cd "$(dirname "$0")/../concordance" && pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
# the benchmark bed: the stratum BB of the concordance strata
BB=$(awk -F'\t' -v cd=$CONC '$1 == "BB" {print (($3 ~ /^\//) ? $3 : cd "/" $3)}' $CONC/config/strata.tsv)

python3 $HERE/lib/truth_preclean.py --truth $TRUTH --ref $REF --bb $BB --out $H/truth_prep --max-in-bb 10
bash $HERE/prepare_callset.sh $H/truth_prep/$(basename $TRUTH .vcf.gz).refclean.vcf.gz HG002 HG002 $H/prep/GIAB_v5.0q/HG002 all chr
