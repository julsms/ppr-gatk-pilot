process PANGENIE_BIALLELIC {
    tag "$sample_name"
    publishDir "${params.outdir}/pangenie", mode: 'symlink'

    // rule convert_back_biallelic_representation of PanGenie's run-from-callset pipeline:
    // the multi-allelic PanGenie genotypes -> one record per ALT allele, using the IDs of the
    // biallelic callset. The PPR paper evaluates this biallelic VCF with truvari (SI 5.1.2).
    container 'ppr/pangenie-callset:4.2.1'

    input:
    tuple val(sample_name), path(genotyping_vcf), path(callset_biallelic), path(callset_tbi)

    output:
    tuple val(sample_name), path("${sample_name}_pangenie.vcf"), emit: vcf

    script:
    def S = '/apps/pangenie/pipelines/run-from-callset/scripts'
    """
    set -o pipefail
    cat ${genotyping_vcf} | python3 ${S}/convert-to-biallelic.py ${callset_biallelic} \\
        > ${sample_name}_pangenie.vcf 2> convert-to-biallelic.log
    """
}
