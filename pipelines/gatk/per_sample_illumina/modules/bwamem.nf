process BWAMEM_SORT {
    tag "$sample_id"
    publishDir "${params.outdir}/aligned_bam/", mode: 'symlink'

    cpus 60
    memory '60 GB'
    time '6.h'

    input:
    tuple val(sample_id), path(read1), path(read2)
    tuple path(reference), path(reference_fai), path(reference_dict), path(reference_amb), path(reference_ann), path(reference_bwt), path(reference_pac), path(reference_sa)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    def read_group = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}_lib\\tPU:${sample_id}_unit"
    def tmp_dir = params.tmp_dir
    // Fixed thread counts: without -K, the bwa batch size depends on the number of threads
    def bwa_cpus = 60
    def view_cpus = 18
    def sort_cpus = 24
    """
    mkdir -p ${tmp_dir}
    echo "Starting BWA-MEM alignment and sorting pipeline..."
    echo "Sample: ${sample_id}"
    echo "Read1: ${read1}"
    echo "Read2: ${read2}"
    echo "Reference: ${reference}"
    echo "CPU allocation - BWA: ${bwa_cpus}, samtools view: ${view_cpus}, samtools sort: ${sort_cpus}"
    
    set -euo pipefail

    bwa mem \\
        -t ${bwa_cpus} \\
        -R "${read_group}" \\
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
    
    echo "Pipeline completed successfully"
    if [[ ! -f "${sample_id}_sorted.bam" ]]; then
        echo "ERROR: Output BAM file was not created"
        exit 1
    fi
    
    echo "Output BAM file size: \$(ls -lh ${sample_id}_sorted.bam)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -E '^Version' | sed 's/Version: //')
        samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_sorted.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: 0.7.17
        samtools: 1.17
    END_VERSIONS
    """
}
