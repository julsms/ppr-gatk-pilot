#!/bin/bash
# Benchmark of one HG002 call set: hap.py v0.3.15 (vcfeval engine) on the whole VCF against GIAB v5.0q, with the benchmark bed
# as -f and the strata in one --stratification file; the same against GIAB CMRG v1.00; then <name>/summary.tsv.
# Usage: bash happy.sh <name> <query VCF>      (outputs in <name>/; resources as in README.md)
set -euo pipefail

NAME=$1
Q=$(readlink -f $2)
HERE=$(cd "$(dirname "$0")" && pwd)
RES=${RES:-$HERE/resources}
V5=$RES/giab_v5.0q                     # HG002_GRCh38_v5.0q_smvar.* and their no_chr copies
SR=$HERE/strata/strata_revision        # region beds and stratification files (not in this release)
REF_CHR=$RES/ref_chr/Homo_sapiens.GRCh38.dna.primary_assembly.fasta     # chr names, as in the PPR pipeline
REF_ENS=$RES/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa            # Ensembl names, as in the GATK pipeline
CMRG_VCF=$SR/sources/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz
CMRG_BED=$SR/sources/HG002_GRCh38_CMRG_smallvar_v1.00.bed
HAPPY=${HAPPY:-hap.py}
RTG=${RTG:-rtg}
THREADS=${THREADS:-16}
export RTG_MEM=12g
mkdir -p $NAME/input $NAME/happy $NAME/cmrg

# chr names (PPR, giraffe): the chr truth, bed, strata and reference; Ensembl names (GATK, BWA-MEM): the no_chr ones
first=$(bcftools view -H $Q | head -n 1 | cut -f1 || true)
if [[ $first == chr* ]]; then
    K=with_chr; REF=$REF_CHR; QC=$Q
    TRUTH=$V5/HG002_GRCh38_v5.0q_smvar.vcf.gz; BED=$V5/HG002_GRCh38_v5.0q_smvar.benchmark.bed
else
    K=no_chr; REF=$REF_ENS; QC=$NAME/input/$(basename $Q .vcf.gz).chr_named.vcf.gz
    TRUTH=$V5/HG002_GRCh38_v5.0q_smvar.no_chr.vcf.gz; BED=$V5/HG002_GRCh38_v5.0q_smvar.benchmark.no_chr.bed
fi

# GIAB v5.0q: BB and the 10 strata from one comparison
$HAPPY $TRUTH $Q -f $BED -r $REF --stratification $SR/strata_ga4gh.$K.tsv \
    --engine=vcfeval --engine-vcfeval-path $RTG --threads $THREADS -o $NAME/happy/$NAME > $NAME/happy/hap.py.log 2>&1

# GIAB CMRG v1.00 has chr names only, so an Ensembl-named VCF is renamed first
if [[ $K == no_chr ]]; then
    awk '{print $1 "\tchr" $1}' $REF.fai | sed 's/\tchrMT$/\tchrM/' > $NAME/input/ensembl_to_chr.map
    bcftools annotate --rename-chrs $NAME/input/ensembl_to_chr.map $Q -Oz -o $QC
    tabix -f -p vcf $QC
fi
$HAPPY $CMRG_VCF $QC -f $CMRG_BED -r $REF_CHR \
    --engine=vcfeval --engine-vcfeval-path $RTG --threads $THREADS -o $NAME/cmrg/$NAME > $NAME/cmrg/hap.py.log 2>&1

# summary.tsv: region, type (SNP, INDEL), filter (ALL, PASS); decompose.py and the figure scripts also look for BENCH_DONE
python3 $HERE/summarize_ga4gh.py $NAME $NAME $SR/strata_ga4gh.with_chr.tsv > $NAME/summary.tsv
touch $NAME/BENCH_DONE
