process VARIANT_ANNOTATOR {
    tag "$sample_name"
    container  'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    input:
    tuple val(sample_name), path(vcf_gz), path(vcf_tbi)
    path reference_fasta
    path reference_fai
    path reference_dict
    path dbsnp
    path dbsnp_tbi
    
    output:
    tuple val(sample_name), path("${sample_name}.annotated.vcf.gz"), emit: vcf
    tuple val(sample_name), path("${sample_name}.annotated.vcf.gz.tbi"), emit: tbi
    
    script:
    """
    # Add all necessary annotations for VQSR, keeping multiallelic variants
    gatk --java-options "-Xms20g -Xmx20g" VariantAnnotator \\
        -R ${reference_fasta} \\
        -V ${vcf_gz} \\
        --dbsnp ${dbsnp} \\
        -O ${sample_name}.annotated.vcf.gz \\
        -A QualByDepth \
        -A RMSMappingQuality \
        -A BaseQualityRankSumTest \
        -A AS_QualByDepth \
        -A MappingQuality \
        -A MappingQualityRankSumTest \
        -A GcContent \
        -A Coverage \
        -A AS_MappingQualityRankSumTest \
        -A ReadPosRankSumTest \
        -A FisherStrand \
        -A StrandOddsRatio \
        -A VariantType \
        -A InbreedingCoeff \
        -A  IndelLength \
        -A AlleleFraction 
    
    # Index with tabix, force overwrite if exists
    tabix -f -p vcf ${sample_name}.annotated.vcf.gz
    """
}
