process FILTER_VARIANTS {
    tag "$sample_name"
    
    input:
    tuple val(sample_name), path(vcf_gz)
    
    output:
    tuple val(sample_name), path("${sample_name}.filtered.vcf.gz"), emit: vcf
    tuple val(sample_name), path("${sample_name}.filtered.vcf.gz.tbi"), emit: tbi
    
    script:
    """
    # Extract only variants (exclude 0/0 and ./.)
    bcftools view -i 'GT!="0/0" && GT!="./."' ${vcf_gz} -O z -o ${sample_name}.filtered.vcf.gz
    
    # Index the filtered VCF with tabix
    tabix -p vcf ${sample_name}.filtered.vcf.gz
    """
}
