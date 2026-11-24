process PANGENIE_INDEX {
    tag "$sample_name"
    
    container "mgibio/pangenie:v4.2.1-bookworm"
    
    input:
    tuple val(sample_name), path(vcf)
    path reference
    
    output:
    tuple val(sample_name), path("pangenie_index*"), emit: index
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    """
    # Run PanGenie index
    PanGenie-index -v ${vcf} -r ${reference} -t ${task.cpus} -o pangenie_index ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangenie: \$(PanGenie --help 2>&1 | grep "version:" | sed 's/version: //' | head -1)
    END_VERSIONS
    """
}
