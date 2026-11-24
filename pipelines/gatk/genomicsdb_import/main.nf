#!/usr/bin/env nextflow
// GenomicsDBImport of the per-interval gVCFs of all samples of one platform: one new store per interval, named as
// pipeline2.nf of the platform reads it (illumina <interval>_genomicsdb, aviti genomicsdb_<interval>).
// The sample sheet has the columns sample and folder; the gVCF of an interval is <folder>/<sample>_<interval>.g.vcf.gz.

nextflow.enable.dsl = 2

params.platform = null
params.samplesheet = null
params.intervals = null
params.genomicsdb_dir = null
params.outdir = 'results'
params.tmp_dir = null

include { GENOMICSDBIMPORT_FROM_GVCFS } from './modules/genomicsdbimport_from_gvcfs.nf'

// The same rule as pipeline2.nf of both platforms
def intervalSafe(String interval) {
    return interval.replaceAll('[^A-Za-z0-9]', '_')
}

def storeName(String interval, String platform) {
    def safe = intervalSafe(interval)
    return platform == 'illumina' ? "${safe}_genomicsdb".toString() : "genomicsdb_${safe}".toString()
}

workflow {
    if (!(params.platform in ['illumina', 'aviti'])) {
        error "Set --platform to illumina or aviti (it selects the store names and the GenomicsDBImport options)"
    }
    ['samplesheet', 'intervals', 'genomicsdb_dir'].each { p ->
        if (!params[p]) {
            error "Parameter --${p} is required"
        }
    }

    // Sample sheet: a header line with the columns sample and folder; other columns are ignored
    def sheetLines = file(params.samplesheet).readLines().findAll { line -> line.trim() && !line.startsWith('#') }
    def header = sheetLines[0].split('\t').collect { h -> h.trim() }
    def iSample = header.indexOf('sample')
    def iFolder = header.indexOf('folder')
    if (iSample < 0 || iFolder < 0) {
        error "The sample sheet needs a header line with the columns sample and folder: ${params.samplesheet}"
    }
    def samples = sheetLines.drop(1).collect { line ->
        def f = line.split('\t')
        [f[iSample].trim(), file(f[iFolder].trim())]
    }
    if (samples.isEmpty()) {
        error "The sample sheet has no sample: ${params.samplesheet}"
    }
    def names = samples.collect { s -> s[0] }
    def duplicates = names.findAll { n -> names.count(n) > 1 }.unique()
    if (duplicates) {
        error "Duplicate sample names in the sample sheet: ${duplicates.join(', ')}"
    }

    def intervals = file(params.intervals).readLines().collect { line -> line.trim() }.findAll { line -> line && !line.startsWith('#') }
    if (intervals.isEmpty()) {
        error "No interval in ${params.intervals}"
    }
    def safeNames = intervals.collect { iv -> intervalSafe(iv) }
    if (safeNames.unique(false).size() != safeNames.size()) {
        error "Two intervals give the same store name in ${params.intervals}"
    }

    // Every gVCF and every index must be there before any import starts
    def missing = []
    samples.each { s ->
        intervals.each { iv ->
            def gvcf = s[1].resolve("${s[0]}_${intervalSafe(iv)}.g.vcf.gz".toString())
            if (!gvcf.exists()) {
                missing << gvcf.toString()
            }
            if (!file("${gvcf}.tbi").exists()) {
                missing << "${gvcf}.tbi".toString()
            }
        }
    }
    if (missing) {
        error "${missing.size()} gVCF or index files are missing, for example:\n  " + missing.take(10).join('\n  ')
    }

    // The import never updates a store: the folder must be new or empty (a resumed run keeps its stores)
    def dbDir = file(params.genomicsdb_dir)
    if (dbDir.exists() && !dbDir.isDirectory()) {
        error "--genomicsdb_dir is not a folder: ${dbDir}"
    }
    if (!workflow.resume && dbDir.exists() && dbDir.list().size() > 0) {
        error "--genomicsdb_dir is not empty: ${dbDir}. Give a new folder (or -resume a run that used it)."
    }
    dbDir.mkdirs()

    log.info "GenomicsDB import from gVCFs: platform ${params.platform}, ${samples.size()} samples, ${intervals.size()} intervals"
    log.info "Stores: ${dbDir}/${storeName(intervals[0], params.platform)} and the others"

    import_ch = channel.fromList(intervals).map { iv ->
        def safe = intervalSafe(iv)
        def gvcfs = samples.collect { s -> s[1].resolve("${s[0]}_${safe}.g.vcf.gz".toString()) }
        def tbis = gvcfs.collect { g -> file("${g}.tbi") }
        tuple(iv, storeName(iv, params.platform), names, gvcfs, tbis)
    }

    GENOMICSDBIMPORT_FROM_GVCFS(import_ch)

    GENOMICSDBIMPORT_FROM_GVCFS.out.check
        .map { _interval, check -> check }
        .collectFile(name: 'import_check.tsv', storeDir: params.outdir, sort: true,
                     seed: "interval\tstore\tsamples_expected\tsamples_in_store\tstatus\n")
        .view { f -> "Import check: ${f}" }
}
