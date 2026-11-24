#!/bin/bash
# The other VQSR levels of Figure 7 (SNP/INDEL 99.0/99.0, 99.7/99.0, 99.9/99.0, 99.9/99.7), from the models of vqsr.sh.
# Run vqsr.sh first, in the same folder.
# Usage: bash vqsr_levels.sh <name> <resource folder>
set -euo pipefail

P=$1
RES=$2

# INDEL 99.7 is not in the INDEL tranche list of vqsr.sh, so ApplyVQSR would cut at the 99.9 tranche.
# The INDEL model is fitted again with 99.7 in the list. The model and the VQSLODs stay the same.
mkdir -p recal_indel_t997
gatk --java-options "-Xms4000m -Xmx12000m" VariantRecalibrator \
    -V $P.sites_only.variant_filtered.vcf.gz -O recal_indel_t997/$P.indels.recal \
    --tranches-file recal_indel_t997/$P.indels.tranches \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.7 -tranche 99.5 -tranche 99.0 -tranche 97.0 \
    -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
    -mode INDEL --max-gaussians 4 \
    -resource:mills,known=false,training=true,truth=true,prior=12 $RES/mills.ensembl.vcf.gz \
    -resource:axiomPoly,known=false,training=true,truth=false,prior=10 $RES/axiomPoly.ensembl.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2 $RES/dbsnp138.ensembl.vcf.gz

# INDEL step, once per INDEL level
for IL in 99.0 99.7; do
    R=.
    [[ $IL == 99.7 ]] && R=recal_indel_t997
    gatk --java-options "-Xms5000m -Xmx6500m" ApplyVQSR \
        -O $P.indel_recalibrated_indel$IL.vcf.gz -V $P.variant_filtered.vcf.gz \
        --recal-file $R/$P.indels.recal --tranches-file $R/$P.indels.tranches \
        --truth-sensitivity-filter-level $IL --create-output-variant-index true -mode INDEL
done

# SNP step for each level pair
for pair in 99.0:99.0 99.7:99.0 99.9:99.0 99.9:99.7; do
    SL=${pair%:*}
    IL=${pair#*:}
    gatk --java-options "-Xms5000m -Xmx6500m" ApplyVQSR \
        -O ${P}_joint_vqsrbp_snp${SL}_indel${IL}.vcf.gz -V $P.indel_recalibrated_indel$IL.vcf.gz \
        --recal-file $P.snps.recal --tranches-file $P.snps.tranches \
        --truth-sensitivity-filter-level $SL --create-output-variant-index true -mode SNP
done
