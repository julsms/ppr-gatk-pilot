process MOSDEPTH {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/qc/mosdepth", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.mosdepth.global.dist.txt"), emit: global_dist
    tuple val(sample_id), path("${sample_id}.mosdepth.region.dist.txt"), emit: region_dist
    tuple val(sample_id), path("${sample_id}.mosdepth.summary.txt"), emit: summary
    tuple val(sample_id), path("${sample_id}.regions.bed.gz"), emit: regions
    tuple val(sample_id), path("${sample_id}.regions.bed.gz.csi"), emit: regions_csi
    tuple val(sample_id), path("${sample_id}.per-base.bed.gz"), optional: true, emit: per_base
    tuple val(sample_id), path("${sample_id}.per-base.bed.gz.csi"), optional: true, emit: per_base_csi
    path "versions.yml", emit: versions

    script:
    """
    echo "=== Mosdepth Coverage Analysis ==="
    echo "Sample: ${sample_id}"
    echo "BAM: ${bam}"
    echo "Window size: ${params.mosdepth_window} bp"
    echo "CPUs: ${task.cpus}"
    echo "==================================="
    
    mosdepth \\
        --threads ${task.cpus} \\
        --by ${params.mosdepth_window} \\
        --no-per-base \\
        --fast-mode \\
        ${sample_id} \\
        ${bam}

    if [[ ! -f "${sample_id}.mosdepth.summary.txt" ]]; then
        echo "ERROR: Mosdepth summary file was not created"
        exit 1
    fi

    echo "Mosdepth completed successfully"
    echo ""
    echo "=== Coverage Summary ==="
    
    # Display summary statistics
    cat ${sample_id}.mosdepth.summary.txt
    
    echo ""
    echo "========================"
    echo ""
    echo "Output files:"
    echo "  - Global distribution: ${sample_id}.mosdepth.global.dist.txt"
    echo "  - Region distribution: ${sample_id}.mosdepth.region.dist.txt"
    echo "  - Summary stats: ${sample_id}.mosdepth.summary.txt"
    
    # Calculate mean coverage from summary
    if [[ -f "${sample_id}.mosdepth.summary.txt" ]]; then
        mean_cov=\$(tail -n 1 ${sample_id}.mosdepth.summary.txt | awk '{print \$4}')
        echo ""
        echo "Mean coverage: \${mean_cov}x"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """
}
