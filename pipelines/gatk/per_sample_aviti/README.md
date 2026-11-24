# Per-sample GATK pipeline, Aviti

One sample per run, from the trimmed read pair. There is no BQSR: Aviti base qualities are well calibrated.

1. BWAMEM_SORT: `bwa mem -t 70` (the task cpus) with the read group, piped into `samtools view -@ 21 -bh` and
   `samtools sort -@ 28` (0.3 and 0.4 of the cpus).
2. Picard MarkDuplicates, then `samtools index`.
3. QC of the duplicate-marked BAM: Picard CollectAlignmentSummaryMetrics, CollectWgsMetrics, CollectGcBiasMetrics,
   QualityScoreDistribution and BamIndexStats; `samtools stats`; `mosdepth --by 500 --no-per-base --fast-mode`.
4. HaplotypeCaller `-ERC GVCF --interval-padding 100 --max-reads-per-alignment-start 0` on each of the 66 intervals,
   then IndexFeatureFile.
5. GenomicsDBImport of the sample on each interval (`--reader-threads 8 --batch-size 50`), into stores that each run adds
   its sample to. The joint runs with HG002 and HG001 do not use these stores (`../genomicsdb_import` imports the gVCFs of
   all samples again); the earlier joint calling of the kinship check read them (`../README.md`).

## Run

    nextflow run main.nf -params-file params.yaml --samplename <sample> \
        --fastq1 <sample>_R1_val_1.fq.gz --fastq2 <sample>_R2_val_2.fq.gz

`params.yaml` gives the paths named in `nextflow.config`: the reference files, `picard`, `genomicsdb_dir` (a new folder)
and `tmp_dir`. Without a profile the tools must be on the `PATH`; with `-profile docker --data_dir <folder>` they run in
the images of `nextflow.config`, and all paths must be under that folder. `markdup_bam_dir` is null as in the cohort runs:
the BAM is then copied to `null/` in the work folder of INDEXMARKDUPBAM.

Outputs in `results/<sample>/`: `bam/` (the duplicate-marked BAM, also the input of BWA-MEM + DeepVariant), `metrics/`
(duplicates), `gvcf/` (66 gVCFs with index), `qc/picard/`, `qc/samtools/`, `qc/mosdepth/`, `pipeline_info/`.

## As run

The 24 cohort genomes, HG002 and HG001 ran these commands: all task scripts of their execution reports match this code.
`main.nf` and the modules were rebuilt from the reports of the cohort runs: the task scripts are exact, and the channel
wiring was inferred from the pipeline graph of the reports. Four cohort genomes ran the Picard QC with 8 cpus and 32 GB
(`-Xmx28g`). The read group has `PL:ILLUMINA`, as in the runs. The cpus and memory in `nextflow.config` enter the
commands (threads and Java heap), and the trace has the 33 columns of the cohort runs.

Tools: Nextflow 25.04.7, BWA-MEM 0.7.17, samtools 1.19.2, Picard 2.26.0, GATK 4.5.0.0, mosdepth 0.3.6. Picard
CollectGcBiasMetrics and QualityScoreDistribution draw their PDF charts with R (`Rscript` on the `PATH`; the GATK image has
R 3.6.2).
