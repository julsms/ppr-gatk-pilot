#!/bin/bash
# One individual of a call set, prepared as the call sets of the concordance step: carried alleles only, split, trimmed and
# left-aligned on the reference, inside the analysis space, no '*', unphased genotypes, and the key table with the strata.
# Usage: bash prepare_callset.sh <VCF> <sample name in the VCF> <name to write> <output prefix> <pass|all> <chr|ensembl>
set -euo pipefail

VCF=$1
NAME=$2
S=$3
OUT=$4
FILTER=$5       # pass: PASS records only; all: every FILTER value
CONTIGS=$6      # ensembl: contig names 1-22, X, Y, renamed to chr names
REF=${REF:-GRCh38_chr.fa}           # GRCh38 with chr names (+ .fai), as in the concordance step
CONC=$(cd "$(dirname "$0")/../concordance" && pwd)
mkdir -p "$(dirname $OUT)"

PASS_OPT=()
[[ $FILTER == pass ]] && PASS_OPT=(-f PASS)
RENAME_OPT=()
[[ $CONTIGS == ensembl ]] && RENAME_OPT=(--rename-chrs $CONC/config/ensembl_to_chr.map)
echo $S > $OUT.sample.txt
awk -F'\t' -v cd=$CONC 'NR > 1 {print $1 "\t" (($3 ~ /^\//) ? $3 : cd "/" $3)}' $CONC/config/strata.tsv > $OUT.strata.tsv

bcftools view -s $NAME -Ov $VCF \
    | bcftools view -a -c 1 -Ov - \
    | bcftools view "${PASS_OPT[@]}" -Ov - \
    | bcftools annotate "${RENAME_OPT[@]}" -x INFO,^FORMAT/GT,FORMAT/GQ -Ov - \
    | bcftools norm -f $REF --check-ref e -m -both --old-rec-tag OLD_REC -Ov - \
    | bcftools view -T $CONC/config/analysis_space.bed -Ov - \
    | bcftools view -i 'GT~"[1-9]"' -Ov - \
    | bcftools view -e 'ALT="*"' -Ov - \
    | bcftools norm -d exact -Ov - \
    | awk -f $CONC/lib/unphase_sort_gt.awk \
    | bcftools reheader -s $OUT.sample.txt - \
    | bcftools view -Oz -o $OUT.vcf.gz -
tabix -f -p vcf $OUT.vcf.gz

# key table: type, SPANDEL, autosome flag and one 0/1 column per stratum (the variant bases overlap the BED)
python3 $CONC/lib/make_keys.py --vcf $OUT.vcf.gz --strata $OUT.strata.tsv --out $OUT.keys.tsv.gz > $OUT.key_counts.tsv
