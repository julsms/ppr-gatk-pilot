# PPR pipeline

Nextflow DSL2 workflow of the personalized pangenome reference method of Sirén et al. (2024), for Nextflow 25.10 or later and
Docker. The runs used `--both`, which also calls variants with vg call and PanGenie; the paper uses only the DeepVariant calls.

## Inputs

- Trimmed paired reads (Trim Galore 0.6.10), for example `<sample>_R1_val_1.fq.gz` and `<sample>_R2_val_2.fq.gz`.
- The HPRC v1.1 Minigraph-Cactus GRCh38 graph and its haplotype index, from
  `s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/`: `hprc-v1.1-mc-grch38.gbz`
  (md5 `39a580dacdb8d7c512b96388f6ffa75a`) and `hprc-v1.1-mc-grch38.hapl` (md5 `5d0a4224c130238cfcd45cf1b28ad79b`).
- GRCh38 with `chr` names and its `.fai` (md5 `582323a46ed6cb4925639ca4ab53706e` and `a368fcfa7065c744241ea5d2fa9ee431`). Its
  source is not recorded. Against the Ensembl primary assembly of the GATK pipeline, `chr` is added to the 194 names and chrY
  is the full chromosome (`chromosome:GRCh38:Y:1:57227415:1`): its PARs hold the chrX PAR sequence, not N.

## Run

    docker build -t ppr/pangenie-callset:4.2.1 docker/pangenie-callset
    docker build -t ppr/vg-samtools:1.52.0-1.17 docker/vg-samtools
    nextflow run main.nf --both --run_deepvariant --reads '<sample>_{R1_val_1,R2_val_2}.fq.gz' --sample_name <sample> \
        --outdir results/<sample> --graph_gbz hprc-v1.1-mc-grch38.gbz --graph_hapl hprc-v1.1-mc-grch38.hapl \
        --reference_fasta <GRCh38 chr>.fasta --reference_fai <GRCh38 chr>.fasta.fai -work-dir work

Outputs in `results/<sample>/`: `small_variants/<sample>.deepvariant.vcf.gz` (the PPR call set),
`small_variants/<sample>_chr_named_rg.bam` (read by `analyses/qc`), `alignment/<sample>.gam`
and `alignment/<sample>_align.err` (read by `analyses/benchmark/decomposition`) and `pipeline_info/` (trace and report).
`alignment/` and `kmc/` hold links into `work/`: copy them before deleting `work/`. The commands of each step are in `modules/`.

Containers (`nextflow.config`): KMC 3.2.4, vg 1.52.0, samtools 1.17 with vg 1.52.0, htslib 1.19.1, vcfbub 0.1.2, vcflib 1.0.15,
PanGenie 4.2.1, DeepVariant 1.9.0 (`google/deepvariant:1.9.0-gpu` on one GPU, `--docker_gpu_opts`; `-profile cpu`: the CPU
image `google/deepvariant:1.9.0`). Up to 70 CPUs and 140 GB per task (`process.resourceLimits` on a smaller machine). The
cohort runs used half the host cores (14-44 threads) and 140 GB per heavy task (`kmc -m132`), 2 threads for vg view, 4 for
the XG index and all cores for DeepVariant (28-88 shards; the CPU image with 80-88 shards for individuals 5, 15 and 22 on
Aviti and 15 and 22 on Illumina); HG002 Aviti as well (24 threads, `kmc -m132`, 48 shards). HG002 Illumina was capped at 64
CPUs (`kmc -t64 -m112`, 40 shards). Illumina 1-8 ran with a single-machine copy of this code (70 threads, 40 shards): the
same small-variant commands apart from `kmc -m120`, the haplotype index and no `--intermediate_results_dir` in run_deepvariant
(PANGENIE_PREP had no bcftools sort). That index (`graph.hapl`, 5,183,927,952 B, md5 `159e986b01160651719a151ff7e1883e`)
was built from the same GBZ with vg 1.52.0; its logs record only `vg index -j graph.dist --snarl-limit 1`. On HG002 the two
indexes sampled the same haplotype sequences (358 paths, identical length and sequence md5 multisets); only the path names
differ, so the sampled GBZ md5 differs.
