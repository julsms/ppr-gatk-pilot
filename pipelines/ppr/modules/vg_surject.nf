process VG_SURJECT {
    tag "$sample_name"
    publishDir "${params.outdir}/small_variants", mode: 'copy'

    // samtools 1.17 + static vg 1.52.0 so surject -> sort runs in one pipe (docker/vg-samtools/Dockerfile).
    container 'ppr/vg-samtools:1.52.0-1.17'

    input:
    tuple val(sample_name), path(gam_file), path(xg_file)

    output:
    tuple val(sample_name), path("${sample_name}_chr_named_rg.bam"), path("${sample_name}_chr_named_rg.bam.bai"), emit: bam
    path "${sample_name}.flagstat.txt", emit: flagstat
    path "${sample_name}_surject.log", emit: log

    script:
    // One pipe, one BAM on disk: surject, chr contig names, read group, sort; then index and flagstat.
    // -i: consecutive GAM records are read pairs.
    def paths        = ((1..22) + ['X', 'Y']).collect { c -> "-p GRCh38#0#chr${c}" }.join(' ')
    def sort_threads = Math.min(task.cpus, 16)
    def sort_mem     = Math.max(1, ((task.memory.toGiga() - 24) / sort_threads).intValue())
    """
    set -o pipefail
    vg surject -x ${xg_file} -b ${gam_file} -t ${task.cpus} --prune-low-cplx -i --max-frag-len 3000 \\
        ${paths} \\
        2> ${sample_name}_surject.log \\
    | samtools reheader -c "sed 's/SN:GRCh38#0#chr/SN:chr/'" - \\
    | samtools addreplacerg -u -r '@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1' -o - - \\
    | samtools sort -@ ${sort_threads} -m ${sort_mem}G -T ${sample_name}.sorttmp -o ${sample_name}_chr_named_rg.bam -

    samtools index -@ ${task.cpus} ${sample_name}_chr_named_rg.bam
    samtools flagstat -@ ${task.cpus} ${sample_name}_chr_named_rg.bam > ${sample_name}.flagstat.txt
    """
}
