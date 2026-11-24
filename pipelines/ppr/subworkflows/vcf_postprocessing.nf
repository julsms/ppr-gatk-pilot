include { FILTER_VARIANTS } from '../modules/filter_variants.nf'
include { VARIANT_ANNOTATOR } from '../modules/variant_annotator.nf'
include { VARIANT_RECALIBRATOR_SNPS } from '../modules/variant_recalibrator_snps.nf'
include { VARIANT_RECALIBRATOR_INDELS } from '../modules/variant_recalibrator_indels.nf'
include { APPLY_VQSR } from '../modules/apply_vqsr.nf'

workflow VCF_POSTPROCESSING {
    take:
    vcf_ch      // tuple val(sample_name), path(vcf_gz)
    
    main:
    // Filter variants (remove 0/0 and ./.)
    FILTER_VARIANTS(vcf_ch)
    
    // Annotate variants
    VARIANT_ANNOTATOR(
        FILTER_VARIANTS.out.vcf.join(FILTER_VARIANTS.out.tbi),
        params.reference_fasta,
        params.reference_fai,
        params.reference_dict,
        params.dbsnp,
        params.dbsnp_tbi
    )
    
    // Variant recalibration for SNPs
    VARIANT_RECALIBRATOR_SNPS(
        VARIANT_ANNOTATOR.out.vcf.join(VARIANT_ANNOTATOR.out.tbi),
        params.reference_fasta,
        params.reference_fai,
        params.reference_dict,
        params.hapmap,
        params.hapmap_tbi,
        params.omni,
        params.omni_tbi,
        params.phase1,
        params.phase1_tbi,
        params.dbsnp,
        params.dbsnp_tbi
    )
    
    // Variant recalibration for INDELs
    VARIANT_RECALIBRATOR_INDELS(
        VARIANT_ANNOTATOR.out.vcf.join(VARIANT_ANNOTATOR.out.tbi),
        params.reference_fasta,
        params.reference_fai,
        params.reference_dict,
        params.mills,
        params.mills_tbi,
        params.dbsnp,
        params.dbsnp_tbi
    )
    
    // Apply VQSR
    // One joined tuple per sample: a process cannot declare the same val twice in strict syntax.
    APPLY_VQSR(
        VARIANT_ANNOTATOR.out.vcf.join(VARIANT_ANNOTATOR.out.tbi)
            .join(VARIANT_RECALIBRATOR_SNPS.out.recal)
            .join(VARIANT_RECALIBRATOR_SNPS.out.tranches)
            .join(VARIANT_RECALIBRATOR_INDELS.out.recal)
            .join(VARIANT_RECALIBRATOR_INDELS.out.tranches),
        params.reference_fasta,
        params.reference_fai,
        params.reference_dict
    )
    
    emit:
    vcf = APPLY_VQSR.out.vcf
    tbi = APPLY_VQSR.out.tbi
}
