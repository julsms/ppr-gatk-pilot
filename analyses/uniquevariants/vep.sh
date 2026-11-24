#!/bin/bash
# VEP 115.2 once on the distinct alleles of all unique sets (cohort and HG002) as a sites-only VCF; the annotation is joined back
# by allele (vep/vep_sites.tsv.gz). Run after extract_unique.sh and hg002_precision.sh.
# Usage: bash vep.sh <VEP cache folder> <reference FASTA> <ClinVar VCF> <dbSNP VCF>
set -euo pipefail

CACHE=$1                            # the folder that holds homo_sapiens/115_GRCh38 (offline cache 115)
FASTA=$2                            # GRCh38 with chr names (+ .fai), the reference of the call sets
CLINVAR=$3                          # clinvar_20260923_GRCh38.chr.vcf.gz
DBSNP=$4                            # dbSNP_b157_GRCh38.chr.vcf.gz
UV_OUT=${UV_OUT:-uniquevariants}
FORK=${FORK:-4}
HERE=$(cd "$(dirname "$0")" && pwd)
V=$UV_OUT/vep
B=$V/batches/1
[[ -e $B/BATCH_DONE ]] && { echo "vep.sh: $B is done; move it away to run VEP again"; exit 1; }
mkdir -p $B

# the sites: every distinct allele of the unique sets, sorted in reference order
python3 $HERE/lib/vep_sites.py sites --roots $UV_OUT/unique,$UV_OUT/hg002/unique --batches $V/batches --fai $FASTA.fai \
    --out-prefix $V/batch
mv $V/batch.sites.vcf $B/sites.vcf
mv $V/batch.all_keys.tsv.gz $V/all_keys.tsv.gz

vep --input_file $B/sites.vcf --output_file $B/vep.vcf --format vcf --vcf --force_overwrite \
    --species homo_sapiens --assembly GRCh38 --fork $FORK --offline --cache --dir_cache $CACHE \
    --fasta $FASTA --everything --canonical --gene_phenotype --variant_class --check_existing --symbol --biotype \
    --custom $CLINVAR,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
    --custom $DBSNP,dbSNP,vcf,exact,0,RS \
    --stats_file $B/stats.html --warning_file $B/warnings.txt
bgzip -f $B/vep.vcf

# one row per allele: the fields read by name from the CSQ header
python3 $HERE/lib/vep_sites.py parse --vep-vcf $B/vep.vcf.gz --out $B/vep_sites.batch.tsv.gz
touch $B/BATCH_DONE                 # merge reads only finished batches
python3 $HERE/lib/vep_sites.py merge --batches $V/batches --all-keys $V/all_keys.tsv.gz --out $V/vep_sites.tsv.gz
