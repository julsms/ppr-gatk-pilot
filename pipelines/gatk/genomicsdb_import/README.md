# GenomicsDB import of all samples of one platform

One GenomicsDBImport per interval, with the gVCFs of all samples in one call. The options are those of the per-sample runs
of the cohort (see `../per_sample_illumina/README.md`, As run):

- Illumina: one `-V` per gVCF, `-L <interval>`, `--tmp-dir`; Java `-XX:+UseParallelGC -XX:ParallelGCThreads=60 -Xmx54g`.
- Aviti: `--sample-name-map`, `-L <interval>`, `--reader-threads 8 --batch-size 50`; Java `-Xmx12g`.

The stores get the names that `pipeline2.nf` of the platform reads: `<interval>_genomicsdb` (Illumina) and
`genomicsdb_<interval>` (Aviti), with each non-alphanumeric character of the interval as `_`. A task checks that each
gVCF holds the sample of its row and that the new store holds all samples of the sheet.

## Run

`samplesheet.tsv` has one row per sample: the sample name of the gVCF header and the folder with its 66 gVCFs
(`<sample>_<interval>.g.vcf.gz` with `.tbi`, as written by the per-sample pipeline). For the paper: the 24 cohort genomes
and HG002, 25 rows. HG001 was imported the same way, with the 24 cohort genomes or alone.

    nextflow run main.nf --platform illumina --samplesheet samplesheet.tsv \
        --genomicsdb_dir <new folder> --outdir import_results --tmp_dir <scratch folder>

Use `--platform aviti` for Aviti (`--tmp_dir` is not used there). `--genomicsdb_dir` must be a new or empty folder.
Without a profile `gatk` must be on the `PATH`; `-profile docker --data_dir <folder>` runs it in the GATK image, with
all paths under that folder.

Outputs: the 66 stores in `--genomicsdb_dir`, and `import_results/import_check.tsv` with one line per interval.

## As run

HG002 and HG001 were imported with this code. In the four execution reports of the imports with the 24 cohort genomes,
all 66 task scripts match it. The cpus and memory of the task only schedule it: the options above are written in the
module.

Tools: Nextflow 25.04.7, GATK 4.5.0.0.
