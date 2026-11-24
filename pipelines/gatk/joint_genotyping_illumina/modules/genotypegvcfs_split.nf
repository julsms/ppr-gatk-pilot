process GENOTYPEGVCFS_SPLIT {
    tag "${interval}"
    publishDir "${params.outdir}/genotypegvcfs_split/", mode: 'symlink'

    cpus 4
    memory '12.GB'
    time '6.h'

    input:
    tuple val(interval), path(genomicsdb)
    tuple path(reference), path(reference_fai), path(reference_dict)

    output:
    tuple val(interval), path("${interval.toString().replaceAll('[^A-Za-z0-9]', '_')}_raw.vcf.gz"), emit: interval_vcf
    path "versions.yml", emit: versions

    script:
    def interval_safe = interval.toString().replaceAll('[^A-Za-z0-9]', '_')
    def tmp_dir = params.tmp_dir
    """
    mkdir -p ${tmp_dir}
    echo "Starting GenotypeGVCFs for interval: ${interval}"
    echo "Input GenomicsDB: ${genomicsdb}"
    echo "Reference: ${reference}"

    set -euo pipefail

    gatk --java-options \\
        "-XX:+UseParallelGC \\
        -XX:ParallelGCThreads=${task.cpus} \\
        -Xmx${task.memory.toGiga()-2}g \\
        -Djava.io.tmpdir=${tmp_dir}" \\
        GenotypeGVCFs \\
        -R ${reference} \\
        -V gendb://${genomicsdb} \\
        -G StandardAnnotation \\
        -G AS_StandardAnnotation \\
        --use-new-qual-calculator \\
        --annotate-with-num-discovered-alleles \\
        -A Coverage \\
        -A QualByDepth \\
        -A FisherStrand \\
        -A StrandOddsRatio \\
        -A MappingQuality \\
        -A MappingQualityRankSumTest \\
        -A ReadPosRankSumTest \\
        -A RMSMappingQuality \\
        -A InbreedingCoeff \\
        -O ${interval_safe}_raw.vcf.gz \\
        --tmp-dir ${tmp_dir}

    echo "GenotypeGVCFs completed successfully for interval: ${interval}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def interval_safe = interval.toString().replaceAll('[^A-Za-z0-9]', '_')
    """
    touch ${interval_safe}_raw.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.4.0.0
    END_VERSIONS
    """
}
