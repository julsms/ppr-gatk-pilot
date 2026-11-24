#!/bin/bash
# The inputs of vep.sh: the Ensembl 115 VEP cache, ClinVar 2026-09-23 with chr contig names, and dbSNP b157 cut to the 25 assembled
# molecules and renamed from their RefSeq accessions with the GRCh38.p14 assembly report.
# Usage: bash vep_resources.sh <output folder>
set -euo pipefail

OUT=$1
ENS=https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache
CV=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38   # the dated file of the run; older files move to the archive folders
SNP=https://ftp.ncbi.nih.gov/snp/latest_release/VCF       # build 157 at the time of the run; older builds are in snp/archive/
ASM=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt
mkdir -p $OUT/cache $OUT/annotations $OUT/downloads
cd $OUT

# the cache, with the BSD sum of the run (as in the Ensembl CHECKSUMS file): cache/homo_sapiens/115_GRCh38
wget -q -O downloads/homo_sapiens_vep_115_GRCh38.tar.gz $ENS/homo_sapiens_vep_115_GRCh38.tar.gz
[[ $(sum downloads/homo_sapiens_vep_115_GRCh38.tar.gz | awk '{print $1+0, $2}') == "1402 24680623" ]]
tar xzf downloads/homo_sapiens_vep_115_GRCh38.tar.gz -C cache/

# ClinVar, with the md5 of the run: 1-22, X and Y to chr1-chr22, chrX and chrY, MT to chrM
wget -q -O downloads/clinvar_20260923.vcf.gz $CV/clinvar_20260923.vcf.gz
wget -q -O downloads/clinvar_20260923.vcf.gz.tbi $CV/clinvar_20260923.vcf.gz.tbi
echo "207bafca2bf725595ad632487f3a9aac  downloads/clinvar_20260923.vcf.gz" | md5sum -c -
{ for c in $(seq 1 22) X Y; do printf '%s\tchr%s\n' $c $c; done; printf 'MT\tchrM\n'; } > downloads/map_clinvar.tsv
bcftools annotate --threads 4 --rename-chrs downloads/map_clinvar.tsv -Oz \
    -o annotations/clinvar_20260923_GRCh38.chr.vcf.gz downloads/clinvar_20260923.vcf.gz
tabix -f -p vcf annotations/clinvar_20260923_GRCh38.chr.vcf.gz

# dbSNP, with the md5 of the run: the 25 assembled molecules, renamed from RefSeq accessions to chr1-chr22, chrX, chrY and chrM
wget -q -O downloads/GCF_000001405.40.gz $SNP/GCF_000001405.40.gz
wget -q -O downloads/GCF_000001405.40.gz.tbi $SNP/GCF_000001405.40.gz.tbi
wget -q -O downloads/GCF_000001405.40_GRCh38.p14_assembly_report.txt $ASM
echo "6a6f313e92a39c337571174dad12cfe1  downloads/GCF_000001405.40.gz" | md5sum -c -
tr -d '\r' < downloads/GCF_000001405.40_GRCh38.p14_assembly_report.txt \
    | awk -F'\t' '!/^#/ && $2=="assembled-molecule" {print $7"\t"$10}' > downloads/map_dbsnp.tsv
touch downloads/GCF_000001405.40.gz.tbi     # so that the index does not look older than the data
bcftools view --threads 4 -r "$(cut -f1 downloads/map_dbsnp.tsv | paste -sd,)" -Ou downloads/GCF_000001405.40.gz \
    | bcftools annotate --threads 4 --rename-chrs downloads/map_dbsnp.tsv -Oz -o annotations/dbSNP_b157_GRCh38.chr.vcf.gz
tabix -f -p vcf annotations/dbSNP_b157_GRCh38.chr.vcf.gz

# md5 of the files that VEP read (a new run gives other values: bcftools writes its command line and the date in the header):
# clinvar_20260923_GRCh38.chr.vcf.gz 2c4f8431dca611945af988e2cb07574c, dbSNP_b157_GRCh38.chr.vcf.gz 153d5d8d1186cd509c01ccfc8f3e9247.
# The assembly report of the run was not kept; the same file, downloaded later, has md5 21f3ac4aa8245a99eb874082051b9dde.
