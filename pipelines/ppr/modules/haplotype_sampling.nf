process HAPLOTYPE_SAMPLING {
    tag "$sample_name"
    publishDir "${params.outdir}/haplotype_sampling", mode: 'copy'

    container 'quay.io/vgteam/vg:v1.52.0'

    input:
    tuple val(sample_name), path(kff_file)
    path full_gbz
    path hapl

    output:
    tuple val(sample_name), path("${sample_name}.gbz"), emit: gbz

    script:
    // PPR method (Siren et al. 2024, SI 2.2): diploid sampling from 32 candidates, reference included.
    """
    vg haplotypes -v 2 -t ${task.cpus} \\
        --include-reference --num-haplotypes 32 --diploid-sampling \\
        -i ${hapl} \\
        -k ${kff_file} \\
        -g ${sample_name}.gbz \\
        ${full_gbz}
    """
}
