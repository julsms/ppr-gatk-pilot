#!/bin/bash
# VQSR of the GATK joint calls as in GATK Best Practices (WARP JointGenotyping v1.7.3, hg38 WGS), GATK 4.5.0.0.
# The call set of the paper: SNP 99.9, INDEL 99.9.
# Usage: bash vqsr.sh <raw joint VCF> <name> <resource folder>
#        EXCESSHET=no bash vqsr.sh ...    the control without the ExcessHet filter (Figure 7c), in its own folder
set -euo pipefail

RAW=$1          # merged_raw_variants.vcf.gz of the joint calling (GenotypeGVCFs)
P=$2            # name of the call set, for example Illumina
RES=$3          # output folder of fetch_resources.sh
SNP_LEVEL=99.9
INDEL_LEVEL=99.9

# 1. ExcessHet hard filter, then a sites-only VCF for the models
FILTER=(--filter-expression "ExcessHet > 54.69" --filter-name ExcessHet)
[[ ${EXCESSHET:-yes} == no ]] && FILTER=()
gatk --java-options "-Xms3000m -Xmx3250m" VariantFiltration "${FILTER[@]}" \
    -O $P.variant_filtered.vcf.gz -V $RAW
gatk --java-options "-Xms3000m -Xmx3250m" MakeSitesOnlyVcf \
    -I $P.variant_filtered.vcf.gz -O $P.sites_only.variant_filtered.vcf.gz

# 2. INDEL model. ApplyVQSR can only cut at a tranche listed here.
gatk --java-options "-Xms4000m -Xmx12000m" VariantRecalibrator \
    -V $P.sites_only.variant_filtered.vcf.gz -O $P.indels.recal --tranches-file $P.indels.tranches \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
    -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
    -mode INDEL --max-gaussians 4 \
    -resource:mills,known=false,training=true,truth=true,prior=12 $RES/mills.ensembl.vcf.gz \
    -resource:axiomPoly,known=false,training=true,truth=false,prior=10 $RES/axiomPoly.ensembl.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2 $RES/dbsnp138.ensembl.vcf.gz

# 3. SNP model. The WGS tranche list of WARP with 99.7 added, so that the 99.7 level of Figure 7 is a true cut.
gatk --java-options "-Xms4000m -Xmx20000m" VariantRecalibrator \
    -V $P.sites_only.variant_filtered.vcf.gz -O $P.snps.recal --tranches-file $P.snps.tranches \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 \
    -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
    -mode SNP --max-gaussians 6 \
    -resource:hapmap,known=false,training=true,truth=true,prior=15 $RES/hapmap.ensembl.vcf.gz \
    -resource:omni,known=false,training=true,truth=true,prior=12 $RES/omni.ensembl.vcf.gz \
    -resource:1000G,known=false,training=true,truth=false,prior=10 $RES/1000G.ensembl.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=7 $RES/dbsnp138.ensembl.vcf.gz

# 4. Apply the INDEL model first, then the SNP model on its output
gatk --java-options "-Xms5000m -Xmx6500m" ApplyVQSR \
    -O $P.indel_recalibrated.vcf.gz -V $P.variant_filtered.vcf.gz \
    --recal-file $P.indels.recal --tranches-file $P.indels.tranches \
    --truth-sensitivity-filter-level $INDEL_LEVEL --create-output-variant-index true -mode INDEL
gatk --java-options "-Xms5000m -Xmx6500m" ApplyVQSR \
    -O ${P}_joint_vqsrbp_snp${SNP_LEVEL}_indel${INDEL_LEVEL}.vcf.gz -V $P.indel_recalibrated.vcf.gz \
    --recal-file $P.snps.recal --tranches-file $P.snps.tranches \
    --truth-sensitivity-filter-level $SNP_LEVEL --create-output-variant-index true -mode SNP
