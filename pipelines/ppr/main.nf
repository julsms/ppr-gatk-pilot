#!/usr/bin/env nextflow

include { KMC_COUNT } from './modules/kmc'
include { HAPLOTYPE_SAMPLING } from './modules/haplotype_sampling'
include { VG_GIRAFFE } from './modules/vg_giraffe'
include { VG_PACK } from './modules/vg_pack'
include { VG_CALL } from './modules/vg_call'
include { VG_DECONSTRUCT } from './modules/vg_deconstruct'
include { VCFBUB_FILTER } from './modules/vcfbub_filter'
include { VCFWAVE } from './modules/vcfwave'
include { PANGENIE_PREP } from './modules/pangenie_prep'
include { PANGENIE_INDEX } from './modules/pangenie_index'
include { PANGENIE_GENOTYPE } from './modules/pangenie_genotype'
include { PANGENIE_BIALLELIC } from './modules/pangenie_biallelic'
include { COMPRESS_INDEX_VCF } from './modules/compress_index'
include { COMPRESS_INDEX_VCF as COMPRESS_INDEX_VCF_PANGENIE } from './modules/compress_index'
include { SMALLVARIANTS_DEEPVARIANT } from './subworkflows/smallvariants_deepvariant'

workflow {
    // Checks live inside the workflow: Nextflow 26 strict syntax rejects top-level statements.
    if (!params.reads)  error "Please provide reads with --reads"
    if (!params.outdir) error "Please provide output directory with --outdir"

    def run_vgcall   = params.vgcall || params.both
    def run_pangenie = params.pangenie || params.both
    if (!run_vgcall && !run_pangenie) {
        error """
        Please specify at least one variant calling method:
        --vgcall    : Run only vg call
        --pangenie  : Run only PanGenie
        --both      : Run both vg call and PanGenie
        """
    }

    // One tuple per read pair. Sample id = file name before the R1/R2 token, unless
    // --sample_name is given (single-sample runs, as published).
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
        .map { id, files -> tuple(params.sample_name ?: id, files) }

    // Graph and reference are staged task inputs; no bind mounts, so any host layout works.
    full_gbz = file(params.graph_gbz, checkIfExists: true)
    hapl     = file(params.graph_hapl, checkIfExists: true)
    ref      = file(params.reference_fasta, checkIfExists: true)
    ref_fai  = file(params.reference_fai, checkIfExists: true)

    KMC_COUNT(reads_ch)
    HAPLOTYPE_SAMPLING(KMC_COUNT.out.kff, full_gbz, hapl)
    gbz_ch = HAPLOTYPE_SAMPLING.out.gbz

    // Every per-sample pairing goes through join(): with several samples in flight, two
    // independent queue channels would pair by arrival order, not by sample. The branches
    // are not serialised; the local executor schedules on declared cpus and memory, so on
    // a 157 GB box two 140 GB tasks still never overlap, while on larger hosts the
    // branches run side by side.
    if (run_vgcall) {
        VG_GIRAFFE(gbz_ch.join(reads_ch))
        // PPR method (Siren et al. 2024, SI Table S4): pack/call on the sampled diploid graph.
        VG_PACK(VG_GIRAFFE.out.gam.join(gbz_ch))
        VG_CALL(VG_PACK.out.pack.join(gbz_ch))
        COMPRESS_INDEX_VCF(VG_CALL.out.vcf)
    }

    // PPR method (SI 5.1.2): deconstruct the sampled diploid graph, vcfbub | vcfwave, then
    // PanGenie's run-from-callset pipeline (merge overlapping alleles into bubbles, genotype,
    // convert back to biallelic for evaluation).
    if (run_pangenie) {
        VG_DECONSTRUCT(gbz_ch)
        VCFBUB_FILTER(VG_DECONSTRUCT.out.vcf)
        VCFWAVE(VCFBUB_FILTER.out.vcf)
        PANGENIE_PREP(VCFWAVE.out.vcf, ref, ref_fai)
        PANGENIE_INDEX(PANGENIE_PREP.out.pangenome_vcf, ref)
        PANGENIE_GENOTYPE(PANGENIE_INDEX.out.index.join(reads_ch))
        PANGENIE_BIALLELIC(PANGENIE_GENOTYPE.out.vcf.join(PANGENIE_PREP.out.callset_biallelic))
        COMPRESS_INDEX_VCF_PANGENIE(PANGENIE_BIALLELIC.out.vcf)
    }

    if (params.run_deepvariant) {
        if (run_vgcall) {
            SMALLVARIANTS_DEEPVARIANT(VG_GIRAFFE.out.gam, gbz_ch, ref, ref_fai)
        }
        else {
            log.warn "DeepVariant requires the vg call branch (GAM files). Skipping DeepVariant."
        }
    }
}
