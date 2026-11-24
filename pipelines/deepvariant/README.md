# BWA-MEM + DeepVariant (DV)

DeepVariant 1.9.0 on the duplicate-marked BWA-MEM alignments of the GATK pipeline, before BQSR. This gives the DV call sets
of the paper: the 23 individuals and HG002, on Illumina and on Aviti. The script was written from the run_info records of the
48 runs, which have the same command and the same three steps.

## Inputs

- `<sample>_markdup.bam` of the GATK pipeline: `pipelines/gatk/per_sample_illumina` run with `--stop_after_markdup true`
  (`results/<sample>/markduplicates/`, with the Picard index `<sample>_markdup.bai`) and `pipelines/gatk/per_sample_aviti`
  (`results/<sample>/bam/`, no index).
- The reference of the GATK pipeline, `Homo_sapiens.GRCh38.dna.primary_assembly.fa` (Ensembl GRCh38, 194 contigs, md5
  `9e6b9465dc708d92bf6d67e9c9fa9389`), with its `.fai` (md5 `d527f3eb6b664020cf4d882b5820056f`).

## Run

    bash deepvariant.sh <sample> <sample>_markdup.bam Homo_sapiens.GRCh38.dna.primary_assembly.fa out/<sample>

1. The index next to the BAM as given (`<BAM>.bai` or the Picard `<sample>_markdup.bai`), else `samtools index -@ 16` (the
   Aviti BAMs). The BAM and the index are mounted as files, so published links into two task folders work.
2. `main_chromosomes.bed`: the 25 contigs 1-22, X, Y and MT of the `.fai`, start 0, end = length, in the order of the `.fai`.
3. The dry run of run_deepvariant (`--dry_run=true`) writes its three commands to `steps.txt`.
4. Each command runs in its own container: make_examples (48 shards) and postprocess_variants without a GPU, call_variants
   on one GPU.

Settings: `--model_type=WGS` for both platforms, `--regions` with the 25 main contigs, `--num_shards=48` (`SHARDS`), the
gVCF with `--output_gvcf`, and the make_examples options of the PPR pipeline:
`--make_examples_extra_args=min_mapping_quality=1,keep_legacy_allele_counter_behavior=true,normalize_reads=true`.

## Outputs

`<sample>.deepvariant.vcf.gz` with `.tbi` (the DV call set), `<sample>.deepvariant.g.vcf.gz` with `.tbi`,
`main_chromosomes.bed`, `steps.txt` and `dv_tmp/` (intermediate files). The sample name in the VCF is the SM of the BAM read
group.

## Tools

DeepVariant 1.9.0, image `google/deepvariant:1.9.0-gpu` at digest
`sha256:e0c8734b8700d945e3ee78d609acb90548f829c874596ffca436af8cf379f87a`; Docker with the NVIDIA Container Toolkit; samtools
1.19.2; awk and sed.
