include { VG_VIEW } from '../modules/vg_view.nf'
include { VG_INDEX_XG } from '../modules/vg_index_xg.nf'
include { VG_SURJECT } from '../modules/vg_surject.nf'
include { DEEPVARIANT } from '../modules/deepvariant.nf'
include { VCF_POSTPROCESSING } from './vcf_postprocessing.nf'

workflow SMALLVARIANTS_DEEPVARIANT {
    take:
    gam_ch      // tuple val(sample_name), path(gam_file)
    gbz_ch      // tuple val(sample_name), path(sampled_gbz)
    ref         // path(reference fasta)
    ref_fai     // path(reference fasta index)

    main:
    // GBZ -> VG -> XG, then surject | reheader | add RG | sort | index in one pipe (one BAM on disk).
    VG_VIEW(gbz_ch)
    VG_INDEX_XG(VG_VIEW.out.vg)
    VG_SURJECT(gam_ch.join(VG_INDEX_XG.out.xg))
    DEEPVARIANT(VG_SURJECT.out.bam, ref, ref_fai)

    vcf_filtered = channel.empty()
    if (params.vcf_postprocessing) {
        VCF_POSTPROCESSING(DEEPVARIANT.out.vcf.map { s, vcf, _tbi -> tuple(s, vcf) })
        vcf_filtered = VCF_POSTPROCESSING.out.vcf
    }

    emit:
    vcf_raw = DEEPVARIANT.out.vcf          // tuple(sample, vcf.gz, tbi)
    vcf_filtered = vcf_filtered
    gvcf = DEEPVARIANT.out.gvcf
    bam = VG_SURJECT.out.bam
}
