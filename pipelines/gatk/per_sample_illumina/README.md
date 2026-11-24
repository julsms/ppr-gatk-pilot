# Per-sample GATK pipeline, Illumina

One sample per run, from the trimmed read pair (Trim Galore `*_val_1.fq.gz` and `*_val_2.fq.gz`):

1. BWAMEM_SORT: `bwa mem -t 60` with the read group, piped into `samtools view -@ 18 -bh` and `samtools sort -@ 24`.
2. Picard MarkDuplicates, then BuildBamIndex.
3. BQSR on each of the 66 intervals: BaseRecalibrator with the known indels and Mills VCFs, GatherBQSRReports,
   ApplyBQSR, then `samtools merge` and `samtools index`. ApplyBQSR writes the reads of each interval, so the reads that
   span one of the 41 internal boundaries are written twice: in HG001, all 34 reads across 50 Mb of chromosome 1 (68 of
   204 records in a 600 bp window), and the same at 100 and 150 Mb; none in two windows inside intervals.
4. HaplotypeCaller `-ERC GVCF --interval-padding 1000 --max-reads-per-alignment-start 0` on each interval, then
   IndexFeatureFile.
5. GenomicsDBImport of the sample on each interval, into stores that each run adds its sample to. The joint runs with
   HG002 and HG001 do not use these stores (`../genomicsdb_import` imports the gVCFs of all samples again); the earlier
   joint calling of the kinship check read them (`../README.md`).

## Run

    nextflow run main.nf -params-file params.yaml --samplename <sample> \
        --fastq1 <sample>_R1_val_1.fq.gz --fastq2 <sample>_R2_val_2.fq.gz

`params.yaml` gives the paths named in `nextflow.config`: the reference files, `picard`, the two known-sites VCFs (Broad
hg38 v0 files with renamed contigs: objects, renaming and md5 in `nextflow.config`), `genomicsdb_persistent_path` (a new
folder) and `tmp_dir`. Without a profile the tools must be on the `PATH`; with `-profile docker --data_dir <folder>` they
run in the images of `nextflow.config`, and all paths must be under that folder. With `--stop_after_markdup true` the run
ends after step 2 and publishes the duplicate-marked BAM; the BAMs of the 23 individuals for BWA-MEM + DeepVariant were
made this way.

Outputs: `results/aligned_bam/` (sorted BAM) and, in `results/<sample>/`, `gvcf/` (66 gVCFs with index), `bqsr/` (BAM
and recalibration table), `markduplicates/` (metrics; the BAM only with `--stop_after_markdup true`), `pipeline_info/`.

## As run

The cpus and memory in `nextflow.config` are those of the runs and enter the commands (Java heap and threads); the bwa
and samtools threads of step 1 are fixed in `modules/bwamem.nf`. HG001 and the duplicate marking of the 23 individuals
ran this code. The per-sample runs of the 24 cohort genomes had the same tool arguments with other Java settings:
BaseRecalibrator and ApplyBQSR `-Xmx6g` (BaseRecalibrator with 5 GC threads in individual 8), GatherBQSRReports
`-Xmx38g`, HaplotypeCaller `-Xmx10g -Xss5m` with 4 threads (20 and 60 in two genomes), IndexFeatureFile `-Xmx1g`,
GenomicsDBImport `-Xmx54g` with 60 GC threads; their gVCFs are those of the joint genotyping. Their execution reports hold
the task scripts of 22 genomes (that of individual 1 lacks 48 GenomicsDBImport tasks). The report of individual 6 holds
only 5 aborted HaplotypeCaller tasks (`-Xmx14g -Xss5m`, 10 threads), and that of individual 11 is empty. HG002 ran
`samtools sort -@ 15 -m 8G` right after bwa, and MarkDuplicates with `-Xms138g -Xmx138g` and 65 GC threads. Tools:
Nextflow 25.04.7, BWA-MEM 0.7.17, samtools 1.19.2, Picard 2.26.0, GATK 4.5.0.0.
