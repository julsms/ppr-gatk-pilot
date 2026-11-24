process VCFWAVE {
    tag "$sample_name"

    // PPR method (Siren et al. 2024, SI 5.1.2): vcfbub -l 0 -a 100000 | vcfwave -I 10000 -t N -n
    // (-n does not exist in any vcflib 1.0.9-1.0.15 release; -I = --inv-min in vcflib).
    container 'quay.io/biocontainers/vcflib:1.0.15--h3fa9d83_1'

    input:
    tuple val(sample_name), path(vcf)

    output:
    tuple val(sample_name), path("${sample_name}_decomposed.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    vcfwave -I 10000 -t ${task.cpus} --quiet ${vcf} > ${sample_name}_decomposed.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfwave: \$(vcfwave --version 2>&1 | head -1 || echo vcflib-1.0.15)
    END_VERSIONS
    """
}
