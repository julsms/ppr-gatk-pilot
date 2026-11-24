process VG_PACK {
    tag "$sample_name"
    publishDir "${params.outdir}/SV_calling", mode: 'symlink'

    container 'quay.io/vgteam/vg:v1.52.0'

    input:
    tuple val(sample_name), path(gam_file), path(gbz_file)

    output:
    tuple val(sample_name), path("${sample_name}.pack"), emit: pack

    script:
    // PPR method (Siren et al. 2024, SI Table S4): pack/call on the SAMPLED diploid graph,
    // not on the full HPRC graph.
    """
    vg pack -x ${gbz_file} -g ${gam_file} -Q 5 -t ${task.cpus} -o ${sample_name}.pack
    """
}
