#!/bin/bash
# The tables of the unique variants: Ti/Tv and het/hom per individual and pair (r27/), the cohort tables (tables/: counts,
# carriers, singletons, consequences, novel alleles, GATK status, RefCall classes) and the HG002 precision (hg002/tables/).
# Usage: bash tables.sh    (after every other step)
set -euo pipefail

export CONC_OUT=${CONC_OUT:-concordance}
export UV_OUT=${UV_OUT:-uniquevariants}
export RELATED=${RELATED:-6}        # one of the two first-degree relatives, left out of the singleton counts
CONC=$(cd "$(dirname "$0")/../concordance" && pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
SAMPLES=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $h["exclude"] != "yes" {print $h["patient"]}' \
    $CONC/config/sample_map.tsv)

# shared, genotype-discordant and unique records of the three pairs of each platform: Ti/Tv and het/hom
for P in illumina aviti; do
    for pr in GATK_PPR GATK_DV DV_PPR; do
        for S in $SAMPLES; do
            python3 $HERE/05_tables.py --titv ${pr}_$P $S
        done
    done
done
python3 $HERE/05_tables.py

# HG002: precision of the unique calls inside the benchmark regions (BB), and check D of the pair classes
for P in illumina aviti; do
    python3 $HERE/lib/truth_status.py tables --usets $HERE/config/hg002_unique_sets.tsv --pairs $HERE/config/hg002_pairs.tsv \
        --platform $P --hg002-out $UV_OUT/hg002 --vep $UV_OUT/vep/vep_sites.tsv.gz --bb-stratum BB
done
