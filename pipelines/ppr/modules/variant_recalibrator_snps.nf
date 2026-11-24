process VARIANT_RECALIBRATOR_SNPS {
    tag "$sample_name"
    container  'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    input:
    tuple val(sample_name), path(vcf_gz), path(vcf_tbi)
    path reference_fasta
    path reference_fai
    path reference_dict
    path hapmap
    path hapmap_tbi
    path omni
    path omni_tbi
    path g1000
    path g1000_tbi
    path dbsnp
    path dbsnp_tbi
    
    output:
    tuple val(sample_name), path("${sample_name}.snps.recal"), emit: recal
    tuple val(sample_name), path("${sample_name}.snps.tranches"), emit: tranches
    tuple val(sample_name), path("${sample_name}.snps.plots.R"), emit: plots
    
    script:
    """
    gatk --java-options "-Xms20g -Xmx20g" VariantRecalibrator \\
        -R ${reference_fasta} \\
        -V ${vcf_gz} \\
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \\
        --resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni} \\
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${g1000} \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \\
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \\
        -mode SNP \\
        --use-allele-specific-annotations \\
        --max-gaussians 6 \\
        -O ${sample_name}.snps.recal \\
        --tranches-file ${sample_name}.snps.tranches 
    """
}
