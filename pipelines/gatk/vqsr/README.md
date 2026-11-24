# VQSR of the GATK joint calls

GATK Best Practices for hg38 whole genomes, as in WARP JointGenotyping v1.7.3 (broadinstitute/warp, tag
`JointGenotyping_v1.7.3`), run with GATK 4.5.0.0. The input is the raw joint VCF of the GATK pipeline
(`merged_raw_variants.vcf.gz`: the 24 cohort genomes and HG002, one file per platform).

The call set of the paper is **SNP 99.9 / INDEL 99.9**.

## Order

1. `bash fetch_resources.sh <reference .dict> resources` downloads the six Broad hg38 resources of `resources.tsv`, checks
   their sizes and md5 values, and renames their contigs to the Ensembl names of the reference (`rename_contigs.awk`).
2. `bash vqsr.sh merged_raw_variants.vcf.gz Illumina resources` (the same for Aviti): ExcessHet filter, sites-only VCF,
   INDEL and SNP models, ApplyVQSR INDEL then SNP at 99.9. Output: `Illumina_joint_vqsrbp_snp99.9_indel99.9.vcf.gz`.
3. `bash extract_sample.sh <joint VCF> <HG002 sample name> HG002.vcf.gz` gives the HG002 calls for the benchmark.
4. Figure 7 only: `bash vqsr_levels.sh Illumina resources` applies the same models at the other four level pairs. Extract
   HG002 from each output (step 3), score it with `analyses/benchmark/happy.sh` and copy its `summary.tsv` into
   `figures/vqsr_levels/data/HG002/` under the name that README gives. `EXCESSHET=no bash vqsr.sh ...` gives the control
   without the ExcessHet filter (Figure 7c, HG001 joint-called with the cohort). Run it in its own folder: it writes the same
   file names as step 2.

## Settings (Supplementary Table S3)

| | INDEL | SNP |
|---|---|---|
| Annotations | FS ReadPosRankSum MQRankSum QD SOR DP | QD MQRankSum ReadPosRankSum FS MQ SOR DP |
| Max Gaussians | 4 | 6 |
| Resources (prior) | Mills (12), Axiom Exome Plus poly (10), dbSNP 138 known (2) | HapMap (15), Omni (12), 1000G phase 1 (10), dbSNP 138 known (7) |
| Tranches | 14 (WARP WGS list) | 13 (WARP WGS list and 99.7) |
| Level of the paper | 99.9 | 99.9 |

The models are not allele-specific (the WGS inputs of WARP set `use_allele_specific_annotations` to false). Sites with
ExcessHet above 54.69 are filtered before VQSR and are not scored.

## Tranches and the filter level

ApplyVQSR does not compute tranches. It reads the tranches file of VariantRecalibrator and cuts at the VQSLOD of the
lowest listed tranche at or above the level. A level that is not in the list gives no error, but a different cut. This is
why 99.7 is in the SNP list, and why `vqsr_levels.sh` fits the INDEL model again with 99.7 in the list. The ApplyVQSR log
names the tranche it used ("Keeping all variants in tranche ...").

## Tools

GATK 4.5.0.0 (the `gatk` launcher needs `python` on the PATH), OpenJDK 17, bcftools 1.20, curl. The run called the GATK
jar with `java -jar`, the options that the `gatk` launcher adds and the heap sizes of the commands. It also set
`-XX:ParallelGCThreads=4`, `-Djava.io.tmpdir`, `--tmp-dir` (and `--TMP_DIR` for MakeSitesOnlyVcf) to a work folder; these do
not change the results and are left out here. All tool arguments are those of the run (checked against its command files and
logs). On chromosome 20 these scripts gave the same record, PASS, ExcessHet and HG002 PASS counts and the same tranche
cut-offs as the chromosome 20 test of the run.
