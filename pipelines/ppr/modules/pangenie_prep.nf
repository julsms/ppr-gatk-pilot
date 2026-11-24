process PANGENIE_PREP {
    tag "$sample_name"
    publishDir "${params.outdir}/pangenie", mode: 'symlink'

    // PanGenie's own run-from-callset pipeline (pipelines/run-from-callset/Snakefile), as
    // required by the PPR paper (SI 5.1.2): validate -> prepare-vcf.py -> add-ids.py ->
    // bcftools norm -m- (biallelic callset) -> merge_vcfs.py (overlapping alleles merged
    // into multi-allelic bubbles = the pangenome VCF PanGenie-index takes).
    // Image = mgibio/pangenie:v4.2.1 + pyfaidx (docker/pangenie-callset/Dockerfile).
    container 'ppr/pangenie-callset:4.2.1'

    input:
    tuple val(sample_name), path(vcf)
    path reference
    path reference_fai

    output:
    tuple val(sample_name), path("${sample_name}_pangenome.vcf"), emit: pangenome_vcf
    tuple val(sample_name), path("${sample_name}_callset-biallelic.vcf.gz"), path("${sample_name}_callset-biallelic.vcf.gz.tbi"), emit: callset_biallelic
    path "versions.yml", emit: versions

    script:
    def S = '/apps/pangenie/pipelines/run-from-callset/scripts'
    """
    set -o pipefail

    # rule validate_vcf: REF alleles must match the reference (-e = error on mismatch)
    bcftools norm --check-ref e --fasta-ref ${reference} ${vcf} > /dev/null 2> validate-vcf.log

    # rule prepare_vcf: drop positions with > 20 % missing alleles
    cat ${vcf} | python3 ${S}/prepare-vcf.py --missing 0.2 2> prepare-vcf.log > input-missing-removed.vcf

    # rule add_ids: one ID per ALT allele. vcfwave realigns alleles and can leave a record behind
    # its predecessor (tabix: "Unsorted positions"), so the callset is sorted before indexing.
    cat input-missing-removed.vcf | python3 ${S}/add-ids.py 2> add-ids.log > callset-unsorted.vcf
    bcftools sort -m 8G -T . -Oz -o callset.vcf.gz callset-unsorted.vcf 2> sort.log
    tabix -p vcf callset.vcf.gz

    # rule normalize: one record per ALT allele
    bcftools norm -m- callset.vcf.gz 2> callset-biallelic.log > ${sample_name}_callset-biallelic.vcf
    bgzip -c ${sample_name}_callset-biallelic.vcf > ${sample_name}_callset-biallelic.vcf.gz
    tabix -p vcf ${sample_name}_callset-biallelic.vcf.gz

    # rule merge_haplotypes: overlapping alleles -> multi-allelic bubbles (the PanGenie panel)
    python3 ${S}/merge_vcfs.py merge -vcf ${sample_name}_callset-biallelic.vcf -r ${reference} -ploidy 2 \\
        2> pangenome.log > ${sample_name}_pangenome.vcf

    rm -f input-missing-removed.vcf callset-unsorted.vcf ${sample_name}_callset-biallelic.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangenie_pipeline: run-from-callset (PanGenie v4.2.1 repo)
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
