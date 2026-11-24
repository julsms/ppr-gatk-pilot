#!/bin/bash
# Benchmark of one HG002 call set on chromosomes 20-22 (Supplementary Figure S3): hap.py v0.3.15 (vcfeval engine) against GIAB v5.0q,
# with the same truth, reference and naming as happy.sh; the v5.0q benchmark bed restricted to chromosomes 20, 21 and 22 is -f, and -l
# limits the comparison to them (a comma-separated list); one comparison, no strata. DeepVariant holds chromosomes 20-22 out of the
# training of its models. STRATA=1 adds the strata file of happy.sh (--stratification, the same 10 strata), so that each region is the
# chromosome 20-22 benchmark region intersected with that stratum, read from the same comparison; outputs then go to <name>/chr20_22_strata/.
# Usage: bash bench_chr20.sh <name> <query VCF>      (the name and VCF of happy.sh; outputs in <name>/chr20_22/; resources as in README.md)
#        STRATA=1 bash bench_chr20.sh <name> <query VCF>    (outputs in <name>/chr20_22_strata/)
set -euo pipefail

NAME=$1
Q=$(readlink -f $2)
HERE=$(cd "$(dirname "$0")" && pwd)
RES=${RES:-$HERE/resources}
V5=$RES/giab_v5.0q
SR=$HERE/strata/strata_revision        # as in happy.sh
REF_CHR=$RES/ref_chr/Homo_sapiens.GRCh38.dna.primary_assembly.fasta
REF_ENS=$RES/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
HAPPY=${HAPPY:-hap.py}
RTG=${RTG:-rtg}
THREADS=${THREADS:-4}
STRATA=${STRATA:-0}
export RTG_MEM=12g
TAG=chr20_22
if [[ $STRATA == 1 ]]; then TAG+=_strata; fi
MARK=${TAG^^}
O=$NAME/$TAG
mkdir -p $O

# chr names (PPR): chr20,chr21,chr22 and the chr truth, bed, strata and reference; Ensembl names (GATK, BWA-MEM): 20,21,22 and the no_chr ones
first=$(bcftools view -H $Q | head -n 1 | cut -f1 || true)
if [[ $first == chr* ]]; then
    K=with_chr; LOC=chr20,chr21,chr22; REF=$REF_CHR
    TRUTH=$V5/HG002_GRCh38_v5.0q_smvar.vcf.gz; BED=$V5/HG002_GRCh38_v5.0q_smvar.benchmark.bed
else
    K=no_chr; LOC=20,21,22; REF=$REF_ENS
    TRUTH=$V5/HG002_GRCh38_v5.0q_smvar.no_chr.vcf.gz; BED=$V5/HG002_GRCh38_v5.0q_smvar.benchmark.no_chr.bed
fi
awk -v l=$LOC 'BEGIN {n = split(l, a, ","); for (i = 1; i <= n; i++) keep[a[i]] = 1} $1 in keep' $BED > $O/benchmark_chr20_22.bed
[[ -s $O/benchmark_chr20_22.bed ]] || { echo "no $LOC lines in $BED"; exit 1; }
SOPT=()
if [[ $STRATA == 1 ]]; then SOPT=(--stratification $SR/strata_ga4gh.$K.tsv); fi

$HAPPY $TRUTH $Q -f $O/benchmark_chr20_22.bed -l $LOC -r $REF ${SOPT[@]+"${SOPT[@]}"} \
    --engine=vcfeval --engine-vcfeval-path $RTG --threads $THREADS -o $O/$NAME.$TAG > $O/hap.py.log 2>&1
touch $O/${MARK}_DONE
