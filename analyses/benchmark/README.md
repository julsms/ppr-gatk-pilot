# HG002 benchmark

Each HG002 call set is compared with the GIAB HG002 v5.0q small-variant benchmark for GRCh38 by `happy.sh`: hap.py v0.3.15 with the
vcfeval engine (rtg-tools 3.12.1) on the whole VCF, the benchmark bed as `-f` and the strata (in `strata/strata_revision/`, not in this release) as
`--stratification`;
then against GIAB CMRG v1.00, which the paper does not report. `summary.tsv` has recall, precision, F1, TP, FN and FP for BB and
each stratum, for SNVs and indels, ALL and PASS.

| folder | what |
|---|---|
| `graph_arms/` | vg giraffe on the GRCh38-only, HPRC default and HPRC d9 graphs, then DeepVariant (Figure 2) |
| `decomposition/` | errors and effects along the five set-ups, alignment metrics |

Run every step from one benchmark folder and call the scripts by their path in the release: each call set gets its folder
there, and `decomposition/` and the figure scripts read them there.

## Call sets (P = Illumina or Aviti)

| name | VCF |
|---|---|
| `P_GATK_HaplotypeCaller_joint_vqsrbp_snp99.9_indel99.9` | HG002 from the GATK joint call set, VQSR 99.9/99.9 (`pipelines/gatk/vqsr/extract_sample.sh`) |
| `P_GATK_Deepvariant_markdup` | BWA-MEM + DeepVariant (BAM of the GATK pipeline, duplicates marked, before BQSR) |
| `P_PPR_Deepvariant` | the PPR pipeline |
| `P_giraffe_GRCh38only_Deepvariant`, `P_giraffe_HPRC_default_Deepvariant`, `P_giraffe_HPRC_d9_Deepvariant` | `graph_arms/giraffe_arm.sh` |

```bash
bash <release>/analyses/benchmark/happy.sh Illumina_PPR_Deepvariant HG002Illumina_SRR14724532.deepvariant.vcf.gz
```

`happy.sh` reads `RES=<folder>` (default `<release>/analyses/benchmark/resources`): `giab_v5.0q/` with the GIAB files
`HG002_GRCh38_v5.0q_smvar.vcf.gz` (and `.tbi`) and `HG002_GRCh38_v5.0q_smvar.benchmark.bed` (used unchanged) and their no_chr
copies; `ref_chr/Homo_sapiens.GRCh38.dna.primary_assembly.fasta` (chr names, PPR) and
`ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa` (Ensembl names, GATK), with `.fai`. `decomposition/` also reads `giab_v4.2.1/`:
`HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` (and `.tbi`) and `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`. A VCF with
Ensembl names (GATK, BWA-MEM) is scored with the no_chr files. The two GATK call sets and `Illumina_giraffe_GRCh38only_Deepvariant`
ran with `THREADS=8`, the others with 16 (same results).

The GIAB files are on `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/`, in `release/AshkenazimTrio/HG002_NA24385_son/`
(each `.tbi` next to its VCF):

| file | md5 |
|---|---|
| `v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz` | `c71acc71069bf7cd7f51eb8fb0c1a1ab`, the file of the run |
| `v5.0q/HG002_GRCh38_v5.0q_smvar.benchmark.bed` | `3366858af85875cbb3c579412ffa4c56`, the file of the run |
| `NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` | `dc750b3807d4af1f7ffec852e9c2f771`; the copy of the run has the same content, compressed differently |
| `NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed` | `97265e922a97c69a0391cf3f92a89b8b`, the bed of the run |

Outputs in `<name>/`: `happy/<name>.*` (hap.py tables, with `roc.Locations.SNP.csv.gz` and `INDEL` for Figure 2a, b),
`cmrg/<name>.*`, `summary.tsv` and the empty `BENCH_DONE` that the later scripts look for. Tools: hap.py v0.3.15 (Python 2 build)
with rtg-tools 3.12.1 and Java 17 (`HAPPY`, `RTG`; `RTG_MEM=12g`), bcftools 1.20, python3.

## Chromosomes 20-22 (Supplementary Figure S3)

The DeepVariant 1.9.0 WGS model of the DV and PPR pipelines was trained on HG002 among other samples, with chromosomes 20-22 held out
of training. `bench_chr20.sh` scores a call set on chromosomes 20, 21 and 22 in one comparison: hap.py as in `happy.sh` (same truth,
reference, naming and resources), with the v5.0q benchmark bed restricted to these chromosomes as `-f` and `-l 20,21,22`
(`chr20,chr21,chr22` for chr names); no strata; `THREADS` 4. `run_chr20_benchmarks.sh` runs it for the six call sets of Figure 1,
one after the other, from the benchmark folder:

```bash
bash <release>/analyses/benchmark/bench_chr20.sh Illumina_PPR_Deepvariant HG002Illumina_SRR14724532.deepvariant.vcf.gz
bash <release>/analyses/benchmark/run_chr20_benchmarks.sh chr20_callsets.tsv    # name <tab> VCF, one line per call set
```

Outputs in `<name>/chr20_22/`: `benchmark_chr20_22.bed` (1,747 intervals, 122,582,881 bp in the run), `<name>.chr20_22.*` (hap.py
tables, with `roc.Locations.{SNP,INDEL}.PASS.csv.gz`) and the empty `CHR20_22_DONE`; `run_chr20_benchmarks.sh` ends with
`CHR20_22_ALL_DONE`.
`figures/benchmark/chr20_table.py` reads them with the whole-genome `summary.tsv`.

By region: with `STRATA=1`, `bench_chr20.sh` adds the strata file of `happy.sh` (`--stratification`, the same 10 strata), so each
region is the chromosome 20-22 benchmark region intersected with that stratum, from the same comparison. Outputs go to
`<name>/chr20_22_strata/` with the marker `CHR20_22_STRATA_DONE`. `STRATA=1 bash run_chr20_benchmarks.sh <list>` runs 12 call sets:
the six of Figure 1 and the three other graph set-ups of Figure 2 on both platforms; it ends with `CHR20_22_STRATA_ALL_DONE`.
`figures/benchmark/chr20_22_strata_table.py` reads them.
