process HAPLOTYPECALLER {
    tag "${sample_id}:${interval}"
    publishDir "${params.outdir}/${sample_id}/gvcf", mode: 'symlink', pattern: "*.g.vcf.gz"

    input:
    tuple val(sample_id), path(bam), path(bai), path(reference), path(reference_fai), path(reference_dict), path(reference_amb), path(reference_ann), path(reference_bwt), path(reference_pac), path(reference_sa), val(interval)

    output:
    tuple val(sample_id), val(interval), path("${sample_id}_${interval.replaceAll('[^A-Za-z0-9]', '_')}.g.vcf.gz"), emit: gvcf
    path "versions.yml", emit: versions

    script:
    def interval_safe = interval.replaceAll('[^A-Za-z0-9]', '_')
    """
    echo "=== HaplotypeCaller ==="
    echo "Sample: ${sample_id}"
    echo "Interval: ${interval}"
    echo "BAM: ${bam}"
    echo "Reference: ${reference}"
    echo "CPUs: ${task.cpus}"
    echo "Memory: ${task.memory.toGiga()}GB"
    echo "======================="
    
    gatk --java-options "-Xmx${task.memory.toGiga() - 2}g -Xss5m" HaplotypeCaller \\
        -R ${reference} \\
        -I ${bam} \\
        -L "${interval}" \\
        --interval-padding 100 \\
        --max-reads-per-alignment-start 0 \\
        --native-pair-hmm-threads ${task.cpus} \\
        -ERC GVCF \\
        -O ${sample_id}_${interval_safe}.g.vcf.gz

    if [[ ! -f "${sample_id}_${interval_safe}.g.vcf.gz" ]]; then
        echo "ERROR: Output GVCF file was not created"
        exit 1
    fi

    echo "HaplotypeCaller completed successfully"
    echo "Output GVCF: ${sample_id}_${interval_safe}.g.vcf.gz"
    echo "GVCF size: \$(ls -lh ${sample_id}_${interval_safe}.g.vcf.gz | awk '{print \$5}')"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """
}
