# Genotype concordance between call sets

Each individual is compared in the nine pairs of `config/pairs.tsv` (A is the baseline side): GATK vs PPR, GATK vs DV and DV vs PPR
on each platform, and each pipeline between Illumina and Aviti; DV is BWA-MEM + DeepVariant. All call sets are prepared the same way.
RTG vcfeval (GA4GH mode) then puts every record of both sides in a class: same allele and genotype (GM), same allele with another
genotype (AM), unique (U) or not assessed (NA). The same records are also matched by exact keys.

## Inputs

- GATK: `<platform>_joint_vqsrbp_snp99.9_indel99.9.vcf.gz` of `pipelines/gatk/vqsr` (one joint VCF per platform). DV:
  `<sample>.deepvariant.vcf.gz` of `pipelines/deepvariant`. PPR: `small_variants/<sample>.deepvariant.vcf.gz` of `pipelines/ppr`.
- `REF`: GRCh38 with chr names (chrMT; the Y PARs are not masked), the reference of the PPR pipeline, with its `.fai`.
- The 11 region BED files (chr names; not in this release); `config/strata.tsv` gives their paths from this folder and their md5.
- `config/sample_map.tsv`: individuals 1-24 (23 is excluded) and, to fill in, the sample name of each one in each call set. This release has no identity check of the pairs (the run had `02_identity_gate.sh`).

## Order (in one output folder named `concordance`; `S=<release>/analyses/concordance`)

`analyses/uniquevariants`, `analyses/qc`, `analyses/kmers`, `figures/qc` and `figures/coverage` read this folder by its name.

1. `bash $S/00_build_sdf.sh $REF ref/GRCh38_chr.sdf`, once.
2. `bash $S/01_prepare_callset.sh <VCF> <sample name in the VCF> <individual> <ensembl|chr> prep/<SET>_<platform>/<individual>`
   for each call set, with `REF` set: `ensembl` for GATK and DV, `chr` for PPR. Output: `.vcf.gz` and `.keys.tsv.gz`.
3. `bash $S/03_vcfeval_pair.sh <pair> <individual>` for each pair: vcfeval, then `classify/<pair>/<individual>.counts.tsv`
   (classes, J_al and J_gt per region and type) and `.records.tsv.gz` (the class of every record).
4. `python3 $S/05_exact_key_level0.py --pair <pair> --patient <individual>` for each pair: `level0/<pair>/<individual>.level0.tsv`
   (exact-key counts) and `.bridge.tsv` (key match against the vcfeval class).
5. `python3 $S/06_summary_stats.py` at the end: `summary/all_counts.tsv` (with the exact-key columns L0_*; the input of
   `figures/concordance`), `genome_counts.tsv`, `pair_means.tsv`, `checks.tsv` and Wilcoxon tables that the paper does not use.

## Notes

- J_al = S / (S + U_A + U_B), J_gt = G / (S + U_A + U_B); S and G are the means of GM + AM and of GM over the two sides.
- vcfeval runs on the whole genome; a record is in a region when its variant bases overlap the BED. CMRG and XPAR are not reported.
- In this layout the `rows_per_pair` check and the identity-gate message of 06 do not apply, and `presplit_PASS_*` stay NA.
- To run a pair again, first remove `vcfeval/<pair>/<individual>*`, `classify/<pair>/<individual>.*` and `level0/<pair>/<individual>.*`:
  04 and 05 read the `.DONE` files (03 writes the one of vcfeval) and report STALE otherwise.

## Tools

bcftools 1.20 (htslib 1.20, tabix 1.21), RTG Tools 3.12.1 on OpenJDK 11.0.1 (`RTG_MEM=12g`, 8 threads; the run started the jar with
`java -Dtalkback=false -Dusage=false -Xmx12g -jar RTG.jar`), Python 3.12.3 with pandas 2.1.4, numpy 1.26.4 and scipy 1.11.4.
