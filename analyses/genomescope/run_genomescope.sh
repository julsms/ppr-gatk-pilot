#!/bin/bash
# GenomeScope 2.0 on the reads of one sample on one platform (main.nf): the reads are subsampled with rasusa to the lower depth
# of the two platforms of the same person, then FastK (k = 21), Histex -G and GenomeScope 2.0.
# Depth: the mosdepth mean coverage of the GATK alignments, mean_coverage of analyses/qc bam_qc_per_sample_<platform>_gatk.tsv
# (an Aviti read set takes the row whose source_label is its number). Target: --target_bases = floor(lower depth x
# 3.1e9), rasusa seed 42.
# Usage: bash run_genomescope.sh <sample> <R1> <R2> <depth of this read set> <depth on the other platform> <output folder>
set -euo pipefail

SAMPLE=$1
R1=$2                       # trimmed reads
R2=$3
COV=$4                      # mean_coverage of this read set (GATK alignments, analyses/qc)
OTHER_COV=$5                # mean_coverage of the same person on the other platform
OUT=$6                      # genomescope_results_illumina or genomescope_results_aviti
GENOME_SIZE=3100000000      # 3.1 Gb
HERE=$(cd "$(dirname "$0")" && pwd)

# Target: the lower depth of the pair x 3.1 Gb; a read set at or below it is not subsampled
TARGET_COV=$(python3 -c "print(min(${COV}, ${OTHER_COV}))")
TARGET_BASES=$(python3 -c "import math; print(int(math.floor(${TARGET_COV} * ${GENOME_SIZE})))")
SUBSAMPLE_ARGS=()
if python3 -c "exit(0 if ${COV} > ${TARGET_COV} else 1)"; then
    SUBSAMPLE_ARGS=(--subsample true --target_bases ${TARGET_BASES})
fi

nextflow run $HERE/main.nf -c $HERE/nextflow.config --sample_id $SAMPLE --fastq_r1 $R1 --fastq_r2 $R2 --output_dir $OUT \
    "${SUBSAMPLE_ARGS[@]}"
