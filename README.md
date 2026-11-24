# ppr-gatk-pilot

Code of the paper "A stepwise pipeline comparison shows improved variant calling with personalized pangenome references across
Illumina and Aviti platforms" (Scientific Reports, submitted).

Three small-variant pipelines were compared: GATK, BWA-MEM with DeepVariant (DV) and the personalized pangenome reference
(PPR) pipeline, on genomes sequenced on Illumina and on Element Aviti, with GIAB HG002 as the benchmark genome and GIAB HG001
as a second check of the GATK VQSR levels.

Each folder holds the commands of one step for one sample (or one cohort step) on one machine, and a README with the order of
the commands, the inputs, the outputs and the tool versions. The commands and their settings are those of the runs behind the
paper. Where a run differed from the code (Java heap or thread settings, for example), the README of the folder says so. Folder
and variable names (for example `strata_revision`, `BENCHMARK_REVISION_DIR` and `kmers_validation/revision`) are those of the runs.

## Samples

24 individuals were sequenced on both platforms and joint-called by GATK (the 24 cohort genomes). Individual 23 was left out of
the analyses for its low Illumina yield (the kinship check keeps all 24), so the paper reports 23 individuals. The analyses
number the individuals 1-24 (23 excluded). Some figure tables number the 23 individuals 1-23, in the text order of the numbers
(1, 10, 11, ..., 9) or in numeric order; each figure README says which.

The sample sheets and maps of `analyses/concordance`, `analyses/genomescope`, `analyses/kmers`, `analyses/performance` and
`figures/coverage` link each individual to its files on each platform.

HG002 is named `HG002_SRR14724532` in the Illumina GATK and DV call sets, `HG002Illumina_SRR14724532` in the Illumina PPR call
set and `HG002Aviti_SRR33551004` in the three Aviti call sets; `analyses/uniquevariants/hg002_precision.sh` expects these names.

## Folders

| Folder | Content |
|---|---|
| `pipelines/gatk/` | BWA-MEM, MarkDuplicates, BQSR, HaplotypeCaller, GenomicsDB import and joint genotyping (Nextflow) |
| `pipelines/gatk/vqsr/` | VQSR as in GATK Best Practices (WARP JointGenotyping v1.7.3), SNP and INDEL 99.9, and the other levels |
| `pipelines/deepvariant/` | DeepVariant 1.9.0 on the duplicate-marked GATK alignments (DV) |
| `pipelines/ppr/` | PPR pipeline (Nextflow): KMC, vg haplotypes, vg giraffe, vg surject, DeepVariant |
| `analyses/benchmark/` | hap.py against GIAB HG002 by region, the graph set-ups and the graph decomposition |
| `analyses/concordance/` | genotype concordance between call sets with RTG vcfeval |
| `analyses/uniquevariants/` | variants found by one pipeline and not by the other |
| `analyses/kmers/` | k-mer support of the calls (Merfin) and the k* copy-number check on HG002 |
| `analyses/genomescope/` | k-mer spectra and GenomeScope 2.0 |
| `analyses/qc/` | read, alignment and call quality; kinship; GATK and PPR coverage |
| `analyses/performance/` | run time, CPU time and memory from the Nextflow traces |
| `figures/` | one folder per topic: the tables and tests that the figures and the text report, and the plots |

## Figures and tables

| In the paper | Analysis | Figure folder |
|---|---|---|
| Figure 1, Supplementary Figures S1-S3 | `analyses/benchmark` (S3: `bench_chr20.sh`) | `figures/benchmark` |
| Figure 2, Supplementary Figures S5-S7, Supplementary Table S1 | `analyses/benchmark` (`graph_arms`, `decomposition`) | `figures/benchmark` |
| Figure 3, Supplementary Figures S8-S12, Supplementary File 4 | `analyses/concordance` | `figures/concordance` |
| Figure 4, Supplementary Figures S4 and S13-S18 | `analyses/uniquevariants` | `figures/uniquevariants` |
| Figure 5, Supplementary Figures S19-S23 | `analyses/kmers` | `figures/merfin` |
| Figure 6, Supplementary Figure S24 | `analyses/qc`, `analyses/kmers`, `analyses/concordance` | `figures/coverage` |
| Figure 7, Supplementary Table S3 | `pipelines/gatk/vqsr` | `figures/vqsr_levels` |
| Supplementary Figures S25-S27 | `analyses/genomescope` | `figures/genomescope` |
| Supplementary Figure S28, Supplementary Note 2 | `analyses/performance` | `figures/performance` |
| Supplementary Table S2, Supplementary Note 1, Supplementary File 2 | `analyses/qc` | `figures/qc` |
| Supplementary File 3 | `analyses/benchmark` (`summary.tsv` and `chr20_22/` of each call set) | `figures/benchmark` |

Supplementary Table S3 is the settings table of `pipelines/gatk/vqsr/README.md`. The Excel files of Supplementary Files 2-4 are
not built here; their values are the tables named in the folder READMEs.

## Order

1. Trim the reads with Trim Galore 0.6.10 (cutadapt 4.4): `trim_galore --paired --retain_unpaired --cores 8 --quality 20
   --length 20 --gzip --fastqc`.
2. Call variants: `pipelines/gatk/` (per sample, GenomicsDB import, joint genotyping, then `vqsr/`), `pipelines/deepvariant/`
   and `pipelines/ppr/`.
3. Run the analyses in this order: `analyses/benchmark`, `analyses/concordance`, `analyses/uniquevariants`, `analyses/kmers`
   and `analyses/qc`, then the others. They read each
   other's outputs in one run folder that holds `concordance/`, `uniquevariants/`, `qc/`, `genomescope/`, `reads/` and
   `kmers_validation/revision/` side by side; the variables in the folder READMEs (for example `CONCORDANCE_DIR`, `UV_DIR` and
   `KV_DIR`) point to them.
4. Make the tables and figures with `figures/<topic>/run_all.sh` (for `figures/qc`, the commands of its README).

## Data

- HG002: Illumina NovaSeq 6000 PCR-free (SRR14724532) and Element Aviti (SRR33551004). HG001: SRR14724533 (Illumina) and
  SRR33551003 (Aviti).
- Truth sets for GRCh38: GIAB HG002 v5.0q; HG002 v4.2.1 in the graph decomposition (the comparison with Sirén et al.); CMRG
  v1.00, which `happy.sh` also runs but the paper does not report; HG001 v4.2.1.
- Regions: the GIAB v3.6 stratifications and the ACMG SF v3.3 and OMIM gene regions. Their BED files are not in this release;
  `analyses/concordance/config/strata.tsv` lists them with their md5.
- References: the Ensembl GRCh38 primary assembly (`Homo_sapiens.GRCh38.dna.primary_assembly.fa`, md5
  `9e6b9465dc708d92bf6d67e9c9fa9389`) for GATK and DV, and GRCh38 with `chr` names (md5 `582323a46ed6cb4925639ca4ab53706e`) for
  PPR, the graph set-ups and the comparisons.
- Graphs: HPRC v1.1 Minigraph-Cactus GRCh38 (`hprc-v1.1-mc-grch38.gbz` and its haplotype index); the d9 graph and the giraffe
  indexes are listed in `analyses/benchmark/graph_arms`.
- The sequencing data of the 24 individuals are available as described in the Data Availability statement of the paper.

## Licence and citation

MIT licence (`LICENSE`). Please cite the paper; `CITATION.cff` has the details.
