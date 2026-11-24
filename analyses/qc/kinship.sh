#!/bin/bash
# Pairwise kinship (KING-robust, PLINK2) of the cohort from the per-sample VCFs of one call set: biallelic autosomal SNPs,
# MAF >= 0.05, at most 5% missing genotypes. KING kinship 0.177-0.354 is first degree, 0.0884-0.177 second degree.
# Usage: bash kinship.sh <name> <per-sample VCF> ...      (for example GATK_illumina vcfs/GATK_illumina/*.GATK.filtered.vcf.gz)
set -euo pipefail

NAME=$1
shift
THREADS=4
mkdir -p $NAME
printf '%s\n' "$@" | sort > $NAME/vcf_list.txt

# one multi-sample VCF; a site that is not in the VCF of a sample is homozygous reference in that sample
bcftools merge --file-list $NAME/vcf_list.txt --missing-to-ref --output-type z --output $NAME/${NAME}_merged.vcf.gz --threads $THREADS
bcftools index --tbi $NAME/${NAME}_merged.vcf.gz

plink2 --vcf $NAME/${NAME}_merged.vcf.gz --autosome --snps-only just-acgt --min-alleles 2 --max-alleles 2 --maf 0.05 --geno 0.05 \
    --make-bed --out $NAME/${NAME}_filtered --threads $THREADS

# <name>_king.kin0: the kinship of every pair
plink2 --bfile $NAME/${NAME}_filtered --make-king-table --out $NAME/${NAME}_king --threads $THREADS
