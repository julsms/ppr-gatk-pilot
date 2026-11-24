process INDEXMARKDUPBAM {
    tag "${sample_id}_markdup"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("${bam}.bai"), emit: bam_bai
    path "versions.yml", emit: versions

    script:
    // with markdup_bam_dir = null (the cohort runs) the copy goes to ./null/ in the work folder
    """
    echo "=== Indexing BAM ==="
    echo "Sample: ${sample_id}"
    echo "BAM file: ${bam}"
    echo "CPUs: ${task.cpus}"
    echo "===================="
    
    samtools index \\
        -@ ${task.cpus} \\
        ${bam}
    
    # Copy BAM and BAI to persistent storage
    mkdir -p ${params.markdup_bam_dir}
    cp ${bam} ${params.markdup_bam_dir}/${bam}
    cp ${bam}.bai ${params.markdup_bam_dir}/${bam}.bai

    if [[ ! -f "${bam}.bai" ]]; then
        echo "ERROR: BAM index file was not created"
        exit 1
    fi

    echo "BAM indexing completed successfully"
    echo "Index file: ${bam}.bai"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
