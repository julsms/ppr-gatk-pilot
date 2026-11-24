# Pipeline-specific (unique) variants

A record of one call set is unique when vcfeval finds no match for it in the other call set of the same individual and platform
(class U of the concordance step). Records with the same allele and another genotype are not unique. There are five sets per
individual and platform: PPR vs GATK, GATK vs PPR, PPR vs DV, DV vs PPR, and PPR vs both. The scripts read the outputs of
`analyses/concordance` (`concordance/prep/`, `concordance/classify/`, the SDF in `concordance/ref/`) and use its config,
`04_classify_count.py` and `lib/`. Run every script from the folder that holds `concordance/`; the outputs go to `uniquevariants/`.

## Order

1. `bash extract_unique.sh <sample> <platform>` for every individual and platform: the unique sets (`unique/`).
2. `bash carriers.sh`: carriers in the full call sets, singletons and one multisample VCF per set (`carriers/`).
3. `bash refcall_classes.sh <sample> <platform> <PPR VCF> <DV VCF>`: RefCall records of the other DeepVariant set (`refcall/`).
4. `bash gatk_status.sh <sample> <platform> <raw joint VCF> <final joint VCF>` for every individual, then
   `bash gatk_status_cohort.sh <platform> <raw joint VCF>`: the GATK status of the PPR-only alleles (`reasons/`).
5. `bash hg002_truth.sh <GIAB v5.0q VCF>` once, then `bash hg002_precision.sh <platform> <GATK joint VCF> <PPR VCF> <DV VCF>`.
6. `bash vep_resources.sh <folder>` once: the VEP cache, and ClinVar and dbSNP with chr names (`<folder>/annotations/`). Then
   `bash vep.sh <folder>/cache <FASTA> <ClinVar VCF> <dbSNP VCF>`: VEP on the distinct alleles of all sets (`vep/`).
7. `bash tables.sh`: Ti/Tv and het/hom (`r27/`), the cohort tables (`tables/`) and the HG002 precision (`hg002/tables/`).
8. In a folder of their own, with `UV_DIR=<uniquevariants>` exported, `python3 novel_profile.py` and `python3 novel_per_person.py`:
   the novel alleles by type, substitution, carriers and region, and per individual (GATK-only; `USET=PPRvsGATK`, `PPRvsDV` and
   `DVvsPPR` for the other sets). `python3 per_person_summary.py GATKvsPPR=per_person.tsv PPRvsGATK=per_person_PPRvsGATK.tsv
   PPRvsDV=per_person_PPRvsDV.tsv DVvsPPR=per_person_DVvsPPR.tsv` gives the per-individual numbers of the Results.
9. `bash xplat_jgt.sh <sample>` for every individual, then `python3 xplat_jgt_summary.py` and `python3 xplat_jgt_novel.py`: the
   genotype concordance (J_gt) between Illumina and Aviti of the unique sets PPR vs GATK, GATK vs PPR, PPR vs DV and DV vs PPR,
   computed as for the full call sets (vcfeval and `04_classify_count.py` of `analyses/concordance`), for all unique variants and
   for the novel ones, with the paired Wilcoxon tests and Holm correction over the 10 regions (`xplat_jgt/tables/`; the
   cross-platform paragraph of the Results).

`prepare_callset.sh` prepares one individual of a call set with the commands of the concordance step (steps 4 and 5).
`figures/uniquevariants` draws Figure 4 and Supplementary Figures S4 and S13-S18 from these tables.

## Settings

- Raw joint VCF: `merged_raw_variants.vcf.gz` of the GATK joint calling; final joint VCF: VQSR SNP 99.9 / INDEL 99.9.
- Singletons: alleles carried by one individual in the full call set of the pipeline. `RELATED` (default 6) is one of the two
  first-degree relatives and is left out.
- VEP: offline cache 115 GRCh38, `--everything`, custom ClinVar 2026-09-23 and dbSNP b157 (exact allele). Novel: absent from
  dbSNP b157, gnomAD genomes and 1000 Genomes, autosomes only.
- vcfeval: `--output-mode=ga4gh --ref-overlap`, 12 GB heap, 8 threads (4 for the HG002 truth runs). HG002 is compared with the
  prepared GIAB v5.0q truth on the whole analysis space; the benchmark regions (stratum BB) come from the key tables of both sides.
- The Python steps read an output only when a `.DONE` file sits next to it. The scripts write them (`carriers.sh` also next to
  the key tables in `concordance/prep/`). Some steps refuse to overwrite their own outputs: remove the outputs of an earlier run
  first. `new` in `reasons/` is the name of the GATK set. `REF` is the reference of the concordance step (`GRCh38_chr.fa`).

## Tools

bcftools 1.20, bgzip and tabix (htslib 1.21), RTG Tools 3.12.1, Ensembl VEP 115.2, Python 3.12 with pandas 2.1.4 and numpy 1.26
(and scipy for `xplat_jgt_summary.py` and `xplat_jgt_novel.py`).
