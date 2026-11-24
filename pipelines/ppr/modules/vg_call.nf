process VG_CALL {
    tag "$sample_name"
    publishDir "${params.outdir}/SV_calling", mode: 'symlink'

    container 'quay.io/vgteam/vg:v1.52.0'

    input:
    tuple val(sample_name), path(pack_file), path(gbz_file)

    output:
    tuple val(sample_name), path("${sample_name}_vgcall.vcf"), emit: vcf
    path "${sample_name}_vgcall.err", emit: log

    script:
    // PPR method (Siren et al. 2024, SI Table S4): call on the SAMPLED diploid graph.
    """
    vg call ${gbz_file} -k ${pack_file} -S GRCh38 -s ${sample_name} -t ${task.cpus} -z \\
        > ${sample_name}_vgcall.vcf 2> ${sample_name}_vgcall.err
    """
}
