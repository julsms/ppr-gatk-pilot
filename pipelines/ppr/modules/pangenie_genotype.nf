process PANGENIE_GENOTYPE {
    tag "$sample_name"
    publishDir "${params.outdir}/pangenie", mode: 'copy'

    container "mgibio/pangenie:v4.2.1-bookworm"

    input:
    tuple val(sample_name), path(index_files), path(reads)

    output:
    tuple val(sample_name), path("${sample_name}_pangenie_genotyping.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    // PanGenie does not read .gz ("INPUT FASTA/Q FILE MUST NOT BE COMPRESSED"). The reads are streamed
    // through a process substitution instead of a decompressed copy on disk (~200 GB Illumina, ~430 GB
    // Aviti). Verified on 1 M Aviti pairs with the real index: genotypes md5-identical.
    """
    PanGenie -f pangenie_index -i <(zcat -f ${reads}) -s ${sample_name} \
        -j ${task.cpus} -t ${task.cpus} -o ${sample_name}_pangenie ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangenie: \$(PanGenie --help 2>&1 | grep "version:" | sed 's/version: //' | head -1)
    END_VERSIONS
    """
}
