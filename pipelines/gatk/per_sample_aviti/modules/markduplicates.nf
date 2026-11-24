process MARKDUPLICATES {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/bam", mode: 'symlink', pattern: "*_markdup.bam"
    publishDir "${params.outdir}/${sample_id}/metrics", mode: 'copy', pattern: "*_markdup.txt"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_markdup.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}_markdup.txt"), emit: metrics
    path "versions.yml", emit: versions

    script:
    def tmp_dir = params.tmp_dir
    """
    mkdir -p ${tmp_dir}
    
    echo "=== Mark Duplicates ==="
    echo "Sample: ${sample_id}"
    echo "Input BAM: ${bam}"
    echo "Memory: ${task.memory.toGiga()}GB"
    echo "CPUs: ${task.cpus}"
    echo "======================="
    
    java \\
        -XX:+UseParallelGC \\
        -XX:ParallelGCThreads=${task.cpus} \\
        -Xmx${task.memory.toGiga() - 2}g \\
        -Djava.io.tmpdir=${tmp_dir} \\
        -jar ${params.picard} MarkDuplicates \\
        I=${bam} \\
        O=${sample_id}_markdup.bam \\
        M=${sample_id}_markdup.txt \\
        VALIDATION_STRINGENCY=LENIENT \\
        ASSUME_SORT_ORDER=coordinate \\
        CREATE_INDEX=false \\
        TMP_DIR=${tmp_dir}
    
    if [[ ! -f "${sample_id}_markdup.bam" ]]; then
        echo "ERROR: Output markdup BAM file was not created"
        exit 1
    fi
    
    echo "MarkDuplicates completed successfully"
    echo "Output BAM size: \$(ls -lh ${sample_id}_markdup.bam | awk '{print \$5}')"
    
    # Display duplicate metrics summary
    if [[ -f "${sample_id}_markdup.txt" ]]; then
        echo ""
        echo "=== Duplicate Metrics Summary ==="
        grep -A 1 "^LIBRARY" ${sample_id}_markdup.txt || true
        echo "================================="
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(java -jar ${params.picard} MarkDuplicates --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """
}
