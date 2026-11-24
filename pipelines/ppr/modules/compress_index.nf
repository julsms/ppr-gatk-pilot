process COMPRESS_INDEX_VCF {
    tag "$sample_name"
    publishDir "${params.outdir}/SV_calling", mode: 'copy'

    container 'quay.io/biocontainers/htslib:1.19.1--h81da01d_1'

    input:
    tuple val(sample_name), path(vcf_file)

    output:
    tuple val(sample_name), path("${vcf_file.baseName}.vcf.gz"),     emit: vcf_gz
    tuple val(sample_name), path("${vcf_file.baseName}.vcf.gz.tbi"), emit: vcf_index

    script:
    // Output keeps the input basename: <sample>_vgcall.vcf.gz, <sample>_pangenie.vcf.gz.
    """
    bgzip -@ ${task.cpus} -c ${vcf_file} > ${vcf_file.baseName}.vcf.gz
    tabix -p vcf ${vcf_file.baseName}.vcf.gz
    """
}
