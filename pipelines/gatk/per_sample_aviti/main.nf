#!/usr/bin/env nextflow
// GATK per-sample pipeline for Aviti reads, without BQSR: BWA-MEM, MarkDuplicates, alignment QC (Picard, samtools stats,
// mosdepth), HaplotypeCaller -ERC GVCF and GenomicsDBImport on each interval of params.intervals. One sample per run.

nextflow.enable.dsl=2

include { BWAMEM_SORT } from './modules/bwamem.nf'
include { MARKDUPLICATES } from './modules/markduplicates.nf'
include { INDEXMARKDUPBAM } from './modules/indexmarkdupbam.nf'
include { PICARD_ALIGNMENT_METRICS } from './modules/picard_alignment_metrics.nf'
include { SAMTOOLS_STATS } from './modules/samtools_stats.nf'
include { MOSDEPTH } from './modules/mosdepth.nf'
include { HAPLOTYPECALLER } from './modules/haplotypecaller.nf'
include { INDEXGVCF } from './modules/indexgvcf.nf'
include { GENOMICSDBIMPORT } from './modules/genomicsdbimport.nf'

workflow {

    println "=" * 80
    println "PIPELINE 1: FASTQ TO GENOMICSDB"
    println "=" * 80
    println "Sample Name    : ${params.samplename}"
    println "FASTQ1         : ${params.fastq1}"
    println "FASTQ2         : ${params.fastq2}"
    println "Output Dir     : ${params.outdir}"
    println "Reference      : ${params.reference}"
    println "Intervals File : ${params.intervals}"
    println "GenomicsDB Dir : ${params.genomicsdb_dir}"
    println "Mosdepth Window: ${params.mosdepth_window} bp"
    println "=" * 80

    reference_ch = Channel.value([
        file(params.reference),
        file(params.reference_fai),
        file(params.reference_dict),
        file(params.reference_amb),
        file(params.reference_ann),
        file(params.reference_bwt),
        file(params.reference_pac),
        file(params.reference_sa)
    ])
    reference = file(params.reference)
    reference_fai = file(params.reference_fai)
    reference_dict = file(params.reference_dict)

    input_ch = Channel.of([params.samplename, file(params.fastq1), file(params.fastq2)])

    intervals_ch = Channel.fromPath(params.intervals)
        .splitText()
        .map { it.trim() }
        .filter { it != "" }

    BWAMEM_SORT(input_ch, reference_ch)
    MARKDUPLICATES(BWAMEM_SORT.out.bam)
    INDEXMARKDUPBAM(MARKDUPLICATES.out.bam)   // [sample_id, bam, bai]

    // QC on the duplicate-marked BAM
    picard_input = INDEXMARKDUPBAM.out.bam_bai.map { sample_id, bam, bai -> [sample_id, bam, bai] }
    PICARD_ALIGNMENT_METRICS(picard_input, reference, reference_fai, reference_dict)

    samtools_input = INDEXMARKDUPBAM.out.bam_bai.map { sample_id, bam, bai -> [sample_id, bam, bai] }
    SAMTOOLS_STATS(samtools_input, reference, reference_fai)

    mosdepth_input = INDEXMARKDUPBAM.out.bam_bai.map { sample_id, bam, bai -> [sample_id, bam, bai] }
    MOSDEPTH(mosdepth_input)

    // No BQSR: HaplotypeCaller reads the duplicate-marked BAM, one task per interval
    haplotypecaller_input = INDEXMARKDUPBAM.out.bam_bai
        .combine(reference_ch)
        .combine(intervals_ch)   // [sample_id, bam, bai, 8 reference files, interval]
    HAPLOTYPECALLER(haplotypecaller_input)

    INDEXGVCF(HAPLOTYPECALLER.out.gvcf)
    GENOMICSDBIMPORT(INDEXGVCF.out.gvcf_tbi)

    GENOMICSDBIMPORT.out.genomicsdb.view { sample_id, interval, db ->
        "GenomicsDB created for ${sample_id}:${interval} -> ${db}"
    }
}

workflow.onComplete {
    println "=" * 80
    println "Pipeline 1 Execution Summary"
    println "=" * 80
    println "Completed at : ${workflow.complete}"
    println "Duration     : ${workflow.duration}"
    println "Success      : ${workflow.success}"
    println "Work Dir     : ${workflow.workDir}"
    println "Output Dir   : ${params.outdir}"
    println "=" * 80
    if (workflow.success) {
        println "Pipeline 1 completed successfully!"
        println "Output Files:"
        println "   - Marked duplicates BAM: ${params.outdir}/${params.samplename}/bam/${params.samplename}_markdup.bam"
        println "   - Duplicate metrics: ${params.outdir}/${params.samplename}/metrics/"
        println "   - QC Reports:"
        println "     - Picard metrics: ${params.outdir}/${params.samplename}/qc/picard/"
        println "     - Samtools stats: ${params.outdir}/${params.samplename}/qc/samtools/"
        println "     - Mosdepth coverage: ${params.outdir}/${params.samplename}/qc/mosdepth/"
        println "   - GVCFs (per interval): ${params.outdir}/${params.samplename}/gvcf/"
        println "   - GenomicsDB workspaces: ${params.genomicsdb_dir}/genomicsdb_<interval>/"
        println "   - Pipeline reports: ${params.outdir}/${params.samplename}/pipeline_info/"
        println "Ready for downstream analysis (genotyping, VQSR, annotation)"
    }
}
