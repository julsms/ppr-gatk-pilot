process APPLY_VQSR {
    tag "$sample_name"
    publishDir "${params.outdir}/small_variants_filtered", mode: 'symlink'

    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'   
 
    input:
    tuple val(sample_name), path(vcf_gz), path(vcf_tbi), path(snps_recal), path(snps_tranches), path(indels_recal), path(indels_tranches)
    path reference_fasta
    path reference_fai
    path reference_dict
    
    output:
    tuple val(sample_name), path("${sample_name}.vqsr_filtered.vcf.gz"), emit: vcf
    tuple val(sample_name), path("${sample_name}.vqsr_filtered.vcf.gz.tbi"), emit: tbi
    
    script:
    """
    # Apply VQSR for SNPs (99.9% sensitivity)
    gatk --java-options "-Xms20g -Xmx20g" ApplyVQSR \\
        -R ${reference_fasta} \\
        -V ${vcf_gz} \\
        --recal-file ${snps_recal} \\
        --tranches-file ${snps_tranches} \\
        -mode SNP \\
        --truth-sensitivity-filter-level 99.9 \\
        --use-allele-specific-annotations \\
        -O ${sample_name}.snps_filtered.vcf.gz
    
    # Apply VQSR for INDELs (99.0% sensitivity)
    gatk --java-options "-Xms20g -Xmx20g" ApplyVQSR \\
        -R ${reference_fasta} \\
        -V ${sample_name}.snps_filtered.vcf.gz \\
        --recal-file ${indels_recal} \\
        --tranches-file ${indels_tranches} \\
        -mode INDEL \\
        --truth-sensitivity-filter-level 99.0 \\
        --use-allele-specific-annotations \\
        -O ${sample_name}.vqsr_filtered.vcf.gz
    
    # Index final VCF with tabix, force overwrite
    tabix -f -p vcf ${sample_name}.vqsr_filtered.vcf.gz
    
    # Clean up intermediate file
    rm ${sample_name}.snps_filtered.vcf.gz*
    """
}
