process INDEXGVCF {
    tag "${sample_id}:${interval}"
    publishDir "${params.outdir}/${sample_id}/gvcf", mode: 'symlink', pattern: "*.tbi"

    input:
    tuple val(sample_id), val(interval), path(gvcf)

    output:
    tuple val(sample_id), val(interval), path(gvcf), path("${gvcf}.tbi"), emit: gvcf_tbi
    path "versions.yml", emit: versions

    script:
    """
    echo "=== Indexing GVCF ==="
    echo "Sample: ${sample_id}"
    echo "Interval: ${interval}"
    echo "GVCF: ${gvcf}"
    echo "====================="
    
    gatk --java-options "-Xmx1g" IndexFeatureFile \\
        -I ${gvcf}

    if [[ ! -f "${gvcf}.tbi" ]]; then
        echo "ERROR: GVCF index file was not created"
        exit 1
    fi

    echo "GVCF indexing completed successfully"
    echo "Index file: ${gvcf}.tbi"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """
}
