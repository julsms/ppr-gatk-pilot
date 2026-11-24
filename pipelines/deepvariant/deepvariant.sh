#!/bin/bash
# DeepVariant 1.9.0 (WGS model) on the duplicate-marked BAM of the GATK pipeline, before BQSR, for one sample.
# Calls on the 25 main contigs of the reference (1-22, X, Y, MT); make_examples with 48 shards, call_variants on one GPU.
# Usage: bash deepvariant.sh <sample> <sample>_markdup.bam <reference .fa> <output folder>
set -euo pipefail

S=$1                    # sample name (the SM of the BAM read group)
BAM=$2                  # <sample>_markdup.bam of the GATK pipeline (a file or the link that Nextflow published)
REF=$(readlink -f $3)   # Homo_sapiens.GRCh38.dna.primary_assembly.fa of the GATK pipeline, with its .fai
OUT=$4
SHARDS=${SHARDS:-48}
IMAGE=google/deepvariant:1.9.0-gpu@sha256:e0c8734b8700d945e3ee78d609acb90548f829c874596ffca436af8cf379f87a
MAIN_CONTIGS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
mkdir -p $OUT/dv_tmp
OUT=$(readlink -f $OUT)
# the BAM and its index are mounted as files, side by side in /input (published links can point into two task folders)
dv() { docker run --rm -u $(id -u):$(id -g) -e HOME=/tmp -v $(readlink -f $BAM):/input/$(basename $BAM):ro \
    -v $(readlink -f $BAI):/input/$(basename $BAI):ro -v $(dirname $REF):/ref:ro -v $OUT:/output -w /output "$@"; }

# 1. The index next to the BAM as given (<BAM>.bai, or the Picard <sample>_markdup.bai); samtools makes one for the Aviti BAMs
BAI=$BAM.bai
[[ -e $BAI ]] || BAI=${BAM%.bam}.bai
[[ -e $BAI ]] || { samtools index -@ 16 $BAM; BAI=$BAM.bai; }

# 2. Regions: the 25 main contigs of the .fai, start 0, end = length, in the order of the .fai
awk -F'\t' -v want="$MAIN_CONTIGS" 'BEGIN { n = split(want, w, " "); for (i = 1; i <= n; i++) ok[w[i]] = 1 }
    ($1 in ok) { print $1 "\t0\t" $2 }' $REF.fai > $OUT/main_chromosomes.bed

# 3. The three commands of run_deepvariant, from its dry run
dv $IMAGE /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=/ref/$(basename $REF) --reads=/input/$(basename $BAM) \
    --output_vcf=/output/$S.deepvariant.vcf.gz --output_gvcf=/output/$S.deepvariant.g.vcf.gz \
    --num_shards=$SHARDS --intermediate_results_dir=/output/dv_tmp --regions=/output/main_chromosomes.bed \
    --make_examples_extra_args=min_mapping_quality=1,keep_legacy_allele_counter_behavior=true,normalize_reads=true \
    --dry_run=true 2>&1 | sed -n 's/^time //p' > $OUT/steps.txt

# 4. make_examples, call_variants (the only step with the GPU) and postprocess_variants, each in its own container
dv $IMAGE bash -c "$(sed -n 1p $OUT/steps.txt)"
dv --gpus 1 $IMAGE bash -c "$(sed -n 2p $OUT/steps.txt)"
dv $IMAGE bash -c "$(sed -n 3p $OUT/steps.txt)"
