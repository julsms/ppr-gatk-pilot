#!/usr/bin/env nextflow
// Joint genotyping of the Illumina GenomicsDB stores: GenotypeGVCFs per interval, then MergeVcfs. Its output
// merged_vcf/merged_raw_variants.vcf.gz is the input of ../vqsr, which makes the GATK call set. The runs also had an
// older VQSR and an annotation step after MergeVcfs; their outputs were not used, and they are left out here.

nextflow.enable.dsl=2

// Parameters for joint calling pipeline
params.genomicsdb_dir = null
params.samplename_joint = null

// Reference
params.reference = null
params.reference_fai = null
params.reference_dict = null

// Intervals file 
params.intervals = null

// Output directory
params.outdir = "results_joint_calling"

// Module includes
include { GENOTYPEGVCFS_SPLIT } from './modules/genotypegvcfs_split.nf'
include { MERGE_INTERVAL_VCFS } from './modules/merge_interval_vcfs.nf'

workflow {
    
    println "=== JOINT CALLING PIPELINE CONFIGURATION ==="
    println "GenomicsDB Directory: ${params.genomicsdb_dir}"
    println "Joint Analysis Name: ${params.samplename_joint}"
    println "Output Dir: ${params.outdir}"
    println "Intervals: ${params.intervals}"
    println "Reference: ${params.reference}"
    println "=============================================="

    
    reference_ch = Channel.value([
        file(params.reference),
        file(params.reference_fai),
        file(params.reference_dict)
    ])

    
    intervals_ch = Channel.fromPath(params.intervals)
        .splitText()
        .map { it.trim() }
        .filter { it != "" }

    
    genomicsdb_input_ch = intervals_ch
        .map { interval ->
            def interval_safe = interval.toString().replaceAll('[^A-Za-z0-9]', '_')
            def genomicsdb_path = file("${params.genomicsdb_dir}/${interval_safe}_genomicsdb")
            
            if (!genomicsdb_path.exists()) {
                error "GenomicsDB not found for interval ${interval}: ${genomicsdb_path}"
            }
            
            [interval, genomicsdb_path]
        }

    genomicsdb_input_ch.count().view { n -> "Found GenomicsDB for ${n} intervals" }

    
    GENOTYPEGVCFS_SPLIT(genomicsdb_input_ch, reference_ch)

    
    interval_vcfs_ch = GENOTYPEGVCFS_SPLIT.out.interval_vcf
        .map { interval, vcf -> vcf }  
        .collect()  

   
    reference_merge_ch = Channel.value([
        file(params.reference),
        file(params.reference_fai), 
        file(params.reference_dict)
    ])

    
    MERGE_INTERVAL_VCFS(interval_vcfs_ch, reference_merge_ch)

    MERGE_INTERVAL_VCFS.out.merged_vcf.view { vcf ->
        "Joint calling done; raw joint VCF for ../vqsr: ${vcf}"
    }
}
