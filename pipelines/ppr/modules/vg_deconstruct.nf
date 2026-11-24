process VG_DECONSTRUCT {
    tag "$sample_name"
    publishDir "${params.outdir}/SV_calling", mode: 'symlink'
    
    container 'quay.io/vgteam/vg:v1.52.0'
    
    input:
    tuple val(sample_name), path(gbz)
    
    output:
    tuple val(sample_name), path("${sample_name}_graph.vcf"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    """
    # Generate VCF from GBZ
    vg deconstruct -P GRCh38 -a ${gbz} -t ${task.cpus} > ${sample_name}_graph.vcf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version | head -n1 | sed 's/vg version //')
    END_VERSIONS
    """
}
