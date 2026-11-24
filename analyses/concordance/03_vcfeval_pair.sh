#!/bin/bash
# RTG vcfeval in GA4GH mode of the two prepared call sets of one individual (A = the baseline side of config/pairs.tsv), then the
# classes and counts per region and type (04_classify_count.py). Run it in the output folder of 01_prepare_callset.sh.
# Usage: bash 03_vcfeval_pair.sh <pair> <individual>
set -euo pipefail

PAIR=$1         # a pair of config/pairs.tsv, for example GATK_PPR_illumina
IND=$2
SDF=${SDF:-ref/GRCh38_chr.sdf}
THREADS=${THREADS:-8}
export RTG_MEM=${RTG_MEM:-12g}
HERE=$(cd "$(dirname "$0")" && pwd)

read A B < <(awk -F'\t' -v p=$PAIR '$1 == p {print $2 "_" $3, $4 "_" $5}' $HERE/config/pairs.tsv)
mkdir -p vcfeval/$PAIR
rtg vcfeval -b prep/$A/$IND.vcf.gz -c prep/$B/$IND.vcf.gz -t $SDF -o vcfeval/$PAIR/$IND \
    --output-mode=ga4gh --ref-overlap -T $THREADS
# 04_classify_count.py reads this file as the sign that vcfeval has finished
touch vcfeval/$PAIR/$IND.DONE
python3 $HERE/04_classify_count.py --pair $PAIR --patient $IND
