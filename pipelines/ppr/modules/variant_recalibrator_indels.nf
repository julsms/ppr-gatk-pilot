process VARIANT_RECALIBRATOR_INDELS {
    tag "$sample_name"
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'   
 
    input:
    tuple val(sample_name), path(vcf_gz), path(vcf_tbi)
    path reference_fasta
    path reference_fai
    path reference_dict
    path mills
    path mills_tbi
    path known_indels
    path known_indels_tbi
    
    output:
    tuple val(sample_name), path("${sample_name}.indels.recal"), emit: recal
    tuple val(sample_name), path("${sample_name}.indels.tranches"), emit: tranches
    tuple val(sample_name), path("${sample_name}.indels.plots.R"), emit: plots
    
    script:
    """
    # Run GATK VariantRecalibrator for INDELs without DP annotation
    gatk --java-options "-Xms20g -Xmx20g" VariantRecalibrator \\
        -R ${reference_fasta} \\
        -V ${vcf_gz} \\
        --resource:mills,known=false,training=true,truth=true,prior=12.0 ${mills} \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${known_indels} \\
        -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \\
        -mode INDEL \\
        --use-allele-specific-annotations \\
        --max-gaussians 4 \\
        -O ${sample_name}.indels.recal \\
        --tranches-file ${sample_name}.indels.tranches \\
        --rscript-file ${sample_name}.indels.plots.R
    """
}
