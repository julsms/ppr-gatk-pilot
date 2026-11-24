# GATK pipeline

GATK best practices for germline short variants on GRCh38, written in Nextflow (DSL2) and run on CPUs. Each sample is
processed on its own, then the samples of one platform are genotyped together.

## Order

1. `per_sample_illumina/` or `per_sample_aviti/`, one run per sample: BWA-MEM, MarkDuplicates, BQSR (Illumina only) and
   HaplotypeCaller `-ERC GVCF` on the 66 intervals of `intervals.list` (1-22, X, Y and MT in blocks of 50 Mb, the last
   block of a chromosome up to 64 Mb). Output: 66 gVCFs per sample, and the duplicate-marked BAM that BWA-MEM +
   DeepVariant also used (Illumina: published only with `--stop_after_markdup true`).
2. `genomicsdb_import/`: one GenomicsDB store per interval with the gVCFs of all samples of one platform. For the paper:
   the 24 cohort genomes and HG002, 25 samples.
3. `joint_genotyping_illumina/` or `joint_genotyping_aviti/`: GenotypeGVCFs on each store, then MergeVcfs. Output: the
   raw joint VCF `merged_vcf/merged_raw_variants.vcf.gz`.
4. `vqsr/`: VQSR of `merged_raw_variants.vcf.gz` as in GATK Best Practices (WARP), SNP 99.9 and INDEL 99.9. This is the
   GATK call set of the paper.

One of the 24 cohort genomes was joint-called but is not among the 23 individuals of the study (excluded for low yield).
HG001 went through the same steps, joint-called with the 24 cohort genomes or alone. The kinship check read per-sample
PASS VCFs made from the annotated joint VCF of an earlier joint calling of the 24 cohort genomes without HG002: GenotypeGVCFs
and MergeVcfs as in step 3, on the GenomicsDB stores built one sample at a time by the per-sample runs, then a VQSR on split
SNP and INDEL files at 99.0 and an annotation step. These last two steps are not in this release.

## Code of each run

| Run | Illumina | Aviti |
|---|---|---|
| per sample, 24 cohort genomes | `per_sample_illumina`, with other Java heap and thread settings (see its README) | `per_sample_aviti` |
| per sample, HG002 | `per_sample_illumina`, with other `samtools sort` and MarkDuplicates settings (see its README) | `per_sample_aviti` |
| per sample, HG001 | `per_sample_illumina` | `per_sample_aviti` |
| raw joint VCF of `vqsr/` | `genomicsdb_import --platform illumina`, then `joint_genotyping_illumina` | `genomicsdb_import --platform aviti`, then `joint_genotyping_aviti` |

Their task scripts (Nextflow execution reports) were compared with this code; the READMEs give the differences.

## Tools

Nextflow 25.04.7, BWA-MEM 0.7.17, samtools 1.19.2, Picard 2.26.0, GATK 4.5.0.0 and mosdepth 0.3.6. The reads were trimmed
before the pipeline with Trim Galore 0.6.10 (cutadapt 4.4): `trim_galore --paired --retain_unpaired --cores 8 --quality 20
--length 20 --gzip --fastqc`, adapter detected by Trim Galore. The
reference is the Ensembl GRCh38 primary assembly with its bwa index: `Homo_sapiens.GRCh38.dna.primary_assembly.fa` (md5
`9e6b9465dc708d92bf6d67e9c9fa9389`) and its `.dict` (md5 `7d1a039228a3cc9093b313ee3e583a07`).
