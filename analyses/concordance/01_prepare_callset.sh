#!/bin/bash
# Prepares one call set of one individual for vcfeval: this sample only, carried ALT alleles, PASS, split and left-aligned on the
# chr-named GRCh38, chr1-22, chrX and chrY without the Y PARs, no '*', genotypes unphased; then its key table (lib/make_keys.py).
# Usage: REF=<FASTA> bash 01_prepare_callset.sh <VCF> <sample name in the VCF> <individual> <ensembl|chr> <output prefix>
set -euo pipefail

IN=$1           # GATK: the joint VCF after VQSR (SNP 99.9 / INDEL 99.9); DV and PPR: the DeepVariant VCF of the individual
LABEL=$2
IND=$3          # sample name of the output, for example 1
CONTIGS=$4      # ensembl for GATK and DV (1, 2, ..., MT), chr for PPR
OUT=$5          # for example prep/GATK_illumina/1
REF=${REF:-GRCh38_chr.fa}       # GRCh38 with chr names and chrMT, Y PARs not masked (the PPR reference), with .fai
HERE=$(cd "$(dirname "$0")" && pwd)

RENAME=()
[[ $CONTIGS == ensembl ]] && RENAME=(--rename-chrs $HERE/config/ensembl_to_chr.map)
mkdir -p $(dirname $OUT)
echo $IND > $OUT.sample.txt

bcftools view -s $LABEL -Ov $IN \
    | bcftools view -a -c 1 -Ov - \
    | bcftools view -f PASS -Ov - \
    | bcftools annotate "${RENAME[@]}" -x INFO,^FORMAT/GT,FORMAT/GQ -Ov - \
    | bcftools norm -f $REF --check-ref e -m -both --old-rec-tag OLD_REC -Ov - \
    | bcftools view -T $HERE/config/analysis_space.bed -Ov - \
    | bcftools view -i 'GT~"[1-9]"' -Ov - \
    | bcftools view -e 'ALT="*"' -Ov - \
    | bcftools norm -d exact -Ov - \
    | awk -f $HERE/lib/unphase_sort_gt.awk \
    | bcftools reheader -s $OUT.sample.txt - \
    | bcftools view -Oz -o $OUT.vcf.gz -
tabix -f -p vcf $OUT.vcf.gz

# key table: type, SPANDEL flag and one 0/1 column per region (the variant bases overlap the BED; relative paths start here)
awk -F'\t' -v cd=$HERE 'NR > 1 {print $1 "\t" (($3 ~ /^\//) ? $3 : cd "/" $3)}' $HERE/config/strata.tsv > $OUT.strata.tsv
python3 $HERE/lib/make_keys.py --vcf $OUT.vcf.gz --strata $OUT.strata.tsv --out $OUT.keys.tsv.gz --bcftools bcftools \
    > $OUT.key_counts.tsv
