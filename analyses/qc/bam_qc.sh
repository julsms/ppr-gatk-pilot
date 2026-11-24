#!/bin/bash
# Coverage and mapping statistics of one BAM (indexed): mosdepth and samtools flagstat.
# Usage: bash bam_qc.sh <BAM> <output prefix>      (for example bam_qc_ppr/illumina_ppr/<sample>)
set -euo pipefail

BAM=$1
OUT=$2
THREADS=${THREADS:-20}
mkdir -p $(dirname $OUT)

# <prefix>.mosdepth.summary.txt and <prefix>.mosdepth.global.dist.txt, no per-base coverage
mosdepth --threads $THREADS --no-per-base --fast-mode $OUT $BAM
samtools flagstat --threads $THREADS $BAM > ${OUT}_flagstat.txt
