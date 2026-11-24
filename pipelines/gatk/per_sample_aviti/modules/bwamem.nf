process BWAMEM_SORT {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(read1), path(read2)
    tuple path(reference), path(reference_fai), path(reference_dict), path(reference_amb), path(reference_ann), path(reference_bwt), path(reference_pac), path(reference_sa)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    def tmp_dir = params.tmp_dir
    def bwa_cpus = task.cpus
    def view_cpus = (task.cpus * 0.3) as int
    def sort_cpus = (task.cpus * 0.4) as int
    """
    mkdir -p ${tmp_dir}
    echo "=== BWA-MEM Alignment Pipeline ==="
    echo "Sample: ${sample_id}"
    echo "Read1: ${read1}"
    echo "Read2: ${read2}"
    echo "Reference: ${reference}"
    echo "CPU allocation - BWA: ${bwa_cpus}, samtools view: ${view_cpus}, samtools sort: ${sort_cpus}"
    echo "================================="
    
    set -euo pipefail

    bwa mem \\
        -t ${bwa_cpus} \\
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}_lib\\tPU:${sample_id}_unit" \\
        "${reference}" \\
        "${read1}" \\
        "${read2}" \\
    | samtools view \\
        -@ ${view_cpus} \\
        -bh \\
        - \\
    | samtools sort \\
        -@ ${sort_cpus} \\
        -T "${tmp_dir}/${sample_id}_tmp" \\
        -o "${sample_id}_sorted.bam" \\
        -
    
    if [[ ! -f "${sample_id}_sorted.bam" ]]; then
        echo "ERROR: Output BAM file was not created"
        exit 1
    fi
    
    echo "BWA-MEM alignment completed successfully"
    echo "Output BAM size: \$(ls -lh ${sample_id}_sorted.bam | awk '{print \$5}')"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -E '^Version' | sed 's/Version: //')
        samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
