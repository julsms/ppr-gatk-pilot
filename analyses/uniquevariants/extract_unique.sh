#!/bin/bash
# The five unique sets of one individual on one platform, from the vcfeval classes of the concordance step: PPR vs GATK,
# GATK vs PPR, PPR vs DV, DV vs PPR, and PPR vs both (unique in the two pairs).
# Usage: bash extract_unique.sh <sample> <platform>
set -euo pipefail

S=$1                                # the individual, as in the concordance sample map
P=$2                                # illumina or aviti
CONC_OUT=${CONC_OUT:-concordance}   # outputs of the concordance step: prep/, classify/
UV_OUT=${UV_OUT:-uniquevariants}
HERE=$(cd "$(dirname "$0")" && pwd)

# extract <set> <side of the pair: A baseline, B calls> <prepared call set of that side> <pair>...
extract() {
    local U=${1}_$P SIDE=$2 SET=$3 rec= cnt= pr
    shift 3
    for pr in "$@"; do
        rec=$rec,$CONC_OUT/classify/${pr}_$P/$S.records.tsv.gz
        cnt=$cnt,$CONC_OUT/classify/${pr}_$P/$S.counts.tsv
    done
    mkdir -p $UV_OUT/unique/$U
    python3 $HERE/lib/extract_unique.py --records ${rec#,} --counts ${cnt#,} --side $SIDE \
        --keys $CONC_OUT/prep/${SET}_$P/$S.keys.tsv.gz --vcf $CONC_OUT/prep/${SET}_$P/$S.vcf.gz \
        --out-prefix $UV_OUT/unique/$U/$S --uset $U --patient $S
    touch $UV_OUT/unique/$U/$S.DONE     # the cohort steps read only the sets with a .DONE file
}

# class U of one side of the pair; PPR vs both: U in both pairs, on the keys of the same PPR file
extract PPRvsGATK B PPR GATK_PPR
extract GATKvsPPR A GATK GATK_PPR
extract PPRvsDV B PPR DV_PPR
extract DVvsPPR A DV DV_PPR
extract PPRvsBOTH B PPR GATK_PPR DV_PPR
