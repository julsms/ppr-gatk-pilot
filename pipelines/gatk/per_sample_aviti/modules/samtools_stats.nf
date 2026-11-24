process SAMTOOLS_STATS {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/qc/samtools", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path(reference)
    path(reference_fai)

    output:
    tuple val(sample_id), path("${sample_id}_samtools_stats.txt"), emit: stats
    path "versions.yml", emit: versions

    script:
    """
    echo "=== Samtools Stats ==="
    echo "Sample: ${sample_id}"
    echo "BAM: ${bam}"
    echo "Reference: ${reference}"
    echo "CPUs: ${task.cpus}"
    echo "======================"
    
    samtools stats \\
        --threads ${task.cpus} \\
        --ref-seq ${reference} \\
        ${bam} \\
        > ${sample_id}_samtools_stats.txt

    if [[ ! -f "${sample_id}_samtools_stats.txt" ]]; then
        echo "ERROR: Samtools stats file was not created"
        exit 1
    fi

    echo "Samtools stats completed successfully"
    echo ""
    echo "=== Key Statistics Summary ==="
    
    # Extract and display key metrics
    grep "^SN" ${sample_id}_samtools_stats.txt | grep -E "raw total sequences|reads mapped|reads mapped and paired|reads duplicated|average length|average quality|insert size average|insert size standard deviation" | \\
        sed 's/^SN\\t/  /' || true
    
    echo "=============================="
    echo "Full stats saved to: ${sample_id}_samtools_stats.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
