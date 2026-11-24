process PICARD_ALIGNMENT_METRICS {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/qc/picard", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path(reference)
    path(reference_fai)
    path(reference_dict)

    output:
    tuple val(sample_id), path("${sample_id}_alignment_metrics.txt"), emit: alignment_metrics
    tuple val(sample_id), path("${sample_id}_wgs_metrics.txt"), emit: wgs_metrics
    tuple val(sample_id), path("${sample_id}_gcbias_metrics.txt"), emit: gcbias_metrics
    tuple val(sample_id), path("${sample_id}_gcbias_summary.txt"), emit: gcbias_summary
    tuple val(sample_id), path("${sample_id}_gcbias_metrics.pdf"), emit: gcbias_pdf
    tuple val(sample_id), path("${sample_id}_quality_score_dist.txt"), emit: quality_score_dist
    tuple val(sample_id), path("${sample_id}_quality_score_dist.pdf"), emit: quality_score_pdf
    tuple val(sample_id), path("${sample_id}_bam_index_stats.txt"), emit: bam_index_stats
    path "versions.yml", emit: versions

    script:
    def tmp_dir = params.tmp_dir
    """
    mkdir -p ${tmp_dir}
    
    echo "=== Comprehensive Picard WGS Metrics ==="
    echo "Sample: ${sample_id}"
    echo "BAM: ${bam}"
    echo "Reference: ${reference}"
    echo "Memory: ${task.memory.toGiga()}GB"
    echo "CPUs: ${task.cpus}"
    echo "========================================"
    
    # 1) CollectAlignmentSummaryMetrics
    echo "[1/6] Running CollectAlignmentSummaryMetrics..."
    java \\
        -XX:+UseParallelGC \\
        -XX:ParallelGCThreads=${task.cpus} \\
        -Xmx${task.memory.toGiga() - 4}g \\
        -Djava.io.tmpdir=${tmp_dir} \\
        -jar ${params.picard} CollectAlignmentSummaryMetrics \\
        R=${reference} \\
        I=${bam} \\
        O=${sample_id}_alignment_metrics.txt \\
        VALIDATION_STRINGENCY=SILENT \\
        TMP_DIR=${tmp_dir}
    
    # 2) CollectWgsMetrics - WGS coverage metrics
    echo "[2/6] Running CollectWgsMetrics..."
    java \\
        -XX:+UseParallelGC \\
        -XX:ParallelGCThreads=${task.cpus} \\
        -Xmx${task.memory.toGiga() - 4}g \\
        -Djava.io.tmpdir=${tmp_dir} \\
        -jar ${params.picard} CollectWgsMetrics \\
        R=${reference} \\
        I=${bam} \\
        O=${sample_id}_wgs_metrics.txt \\
        VALIDATION_STRINGENCY=SILENT \\
        TMP_DIR=${tmp_dir}
    
    # 3) CollectGcBiasMetrics - GC bias + PDF chart
    echo "[3/6] Running CollectGcBiasMetrics..."
    java \\
        -XX:+UseParallelGC \\
        -XX:ParallelGCThreads=${task.cpus} \\
        -Xmx${task.memory.toGiga() - 4}g \\
        -Djava.io.tmpdir=${tmp_dir} \\
        -jar ${params.picard} CollectGcBiasMetrics \\
        R=${reference} \\
        I=${bam} \\
        O=${sample_id}_gcbias_metrics.txt \\
        S=${sample_id}_gcbias_summary.txt \\
        CHART=${sample_id}_gcbias_metrics.pdf \\
        VALIDATION_STRINGENCY=SILENT \\
        TMP_DIR=${tmp_dir}
    
    # 4) QualityScoreDistribution - Quality distribution + PDF
    echo "[4/6] Running QualityScoreDistribution..."
    java \\
        -XX:+UseParallelGC \\
        -XX:ParallelGCThreads=${task.cpus} \\
        -Xmx${task.memory.toGiga() - 4}g \\
        -Djava.io.tmpdir=${tmp_dir} \\
        -jar ${params.picard} QualityScoreDistribution \\
        INPUT=${bam} \\
        OUTPUT=${sample_id}_quality_score_dist.txt \\
        CHART=${sample_id}_quality_score_dist.pdf \\
        VALIDATION_STRINGENCY=SILENT
    
    # 5) BamIndexStats - Per-contig counts
    echo "[5/6] Running BamIndexStats..."
    java \\
        -XX:+UseParallelGC \\
        -XX:ParallelGCThreads=${task.cpus} \\
        -Xmx${task.memory.toGiga() - 4}g \\
        -Djava.io.tmpdir=${tmp_dir} \\
        -jar ${params.picard} BamIndexStats \\
        INPUT=${bam} \\
        > ${sample_id}_bam_index_stats.txt
    
    echo "[6/6] Validating outputs..."
    
    # Validate all outputs
    for file in ${sample_id}_alignment_metrics.txt \\
                ${sample_id}_wgs_metrics.txt \\
                ${sample_id}_gcbias_metrics.txt \\
                ${sample_id}_quality_score_dist.txt \\
                ${sample_id}_bam_index_stats.txt; do
        if [[ ! -f "\$file" ]]; then
            echo "ERROR: \$file was not created"
            exit 1
        fi
    done

    echo ""
    echo "All Picard metrics completed successfully"
    echo ""
    echo "=== Key Metrics Summary ==="
    
    # Display alignment summary
    if grep -q "^PAIR" ${sample_id}_alignment_metrics.txt; then
        echo "Alignment Metrics:"
        grep -A 1 "^PAIR" ${sample_id}_alignment_metrics.txt | tail -n 1 | \\
            awk '{printf "  Total Reads: %s\\n  Aligned: %s (%s%%)\\n  Mean Read Length: %s\\n", \$2, \$7, \$8*100, \$17}'
    fi
    
    # Display WGS coverage summary
    if grep -q "^GENOME_TERRITORY" ${sample_id}_wgs_metrics.txt; then
        echo ""
        echo "WGS Coverage Metrics:"
        grep -A 1 "^GENOME_TERRITORY" ${sample_id}_wgs_metrics.txt | tail -n 1 | \\
            awk '{printf "  Mean Coverage: %.2fx\\n  Median Coverage: %sx\\n  PCT 10x: %s%%\\n  PCT 20x: %s%%\\n  PCT 30x: %s%%\\n", \$2, \$3, \$15*100, \$16*100, \$17*100}'
    fi
    
    echo "============================="
    echo ""
    echo "Output files:"
    echo "  - Alignment metrics: ${sample_id}_alignment_metrics.txt"
    echo "  - WGS coverage: ${sample_id}_wgs_metrics.txt"
    echo "  - GC bias: ${sample_id}_gcbias_metrics.txt + PDF"
    echo "  - Quality dist: ${sample_id}_quality_score_dist.txt + PDF"
    echo "  - Index stats: ${sample_id}_bam_index_stats.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(java -jar ${params.picard} CollectAlignmentSummaryMetrics --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """
}
