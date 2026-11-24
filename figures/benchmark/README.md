# HG002 benchmark figures

Figure 1 and Supplementary Figures S1 to S3 (GATK, BWA-MEM + DeepVariant and PPR; S3 on chromosomes 20-22), Figure 2,
Supplementary Figures S5, S6 and S7 and Supplementary Table S1 (the graph decomposition). The scripts read the outputs of
`analyses/benchmark/` in the benchmark folder (`BENCHMARK_REVISION_DIR`), write their tables next to the scripts and the figures into
`figures/`.

`bash <release>/figures/benchmark/run_all.sh <benchmark folder>` runs them in the order of the run:

```bash
python3 benchmark_table.py                        # benchmark_data.tsv
python3 figure1_panels.py                         # Figure 1a, 1b as single panels
python3 figure1_combine.py                        # Figure 1
python3 supp_metrics_by_region.py                 # Supplementary Figure S1 (PASS)
python3 supp_metrics_by_region.py --filter ALL    # Supplementary Figure S2 (GATK before VQSR)
python3 chr20_table.py                            # chr20_22_benchmark.tsv, chr20_22_roc.tsv
python3 supp_chr20_pr.py                          # Supplementary Figure S3 (chromosomes 20-22)
python3 chr20_compare.py                          # the chromosome 20-22 numbers of the Results
python3 graph_arms_table.py                       # graph_arms_data.tsv, graph_arms_roc.tsv, graph_arms_mapping.tsv
python3 figure2_panels.py                         # Figure 2a to 2e as single panels
python3 figure2_combine.py                        # Figure 2
python3 supp_graph_arms_f1.py                     # Supplementary Figures S5 (Illumina) and S6 (Aviti)
python3 supp_graph_arms_alignment.py              # Supplementary Figure S7
python3 supp_table_graph_setups_errors.py         # Supplementary Table S1 (graph_setups_errors.tsv)
```

| figure or table | file |
|---|---|
| Figure 1, Figure 2 | `figures/Figure1.{pdf,png,tiff}`, `figures/Figure2.{pdf,png,tiff}` |
| S1, S2 | `figures/SuppFig_metrics_by_region_PASS.*`, `figures/SuppFig_metrics_by_region_ALL.*` |
| S3 | `figures/SuppFig_chr20_22_PR.{pdf,png,tiff}` |
| S5, S6 | `figures/SuppFig_graph_arms_F1_Illumina.*`, `figures/SuppFig_graph_arms_F1_Aviti.*` |
| S7 | `figures/SuppFig_graph_arms_alignment.*` |
| Table S1 | `graph_setups_errors.tsv` (errors in BB and the change against BWA-MEM; the script also prints PPR against HPRC default) |

Inputs: `<call set>/summary.tsv` and `BENCH_DONE` of the call sets listed in `analyses/benchmark/README.md`,
`<call set>/happy/<call set>.roc.Locations.{SNP,INDEL}.csv.gz` (Figure 2a, b), `mapping/vg_stats_a.txt` of the giraffe arms and of
`<Platform>_PPR_mapping_stats/`, `comparisons/Mapping_metrics_graph_arms_<Platform>/mapping_metrics.tsv`
(`analyses/benchmark/decomposition/`), and the BAM metrics of BWA-MEM in `comparisons/Mapping_metrics_graph_arms_<Platform>/bam/`
(all contigs) and `bam_chr1_22XY/` (chr1-22, X, Y), made by `analyses/benchmark/decomposition/bam_metrics.sh`. Figure 2e shows
reads mapped and reads with MAPQ 60 in the BAM files that DeepVariant read, on chr1-22, X and Y, for all five set-ups (the giraffe
BAMs after surjection, BWA-MEM counted on the same chromosomes); Supplementary Figure S7 shows all eight metrics, and in its GAM
aligned and GAM MAPQ 60 panels the BWA-MEM row shows its BAM on all contigs. Supplementary Figure S3 reads
`<call set>/chr20_22/` of the six call sets of Figure 1 (`analyses/benchmark/bench_chr20.sh`). `graph_setups.py` holds the names,
order and colours of the five set-ups.
`graph_arms_table.py` stops if F1 does not match precision and recall, if the ROC point at the lowest QUAL threshold is not the PASS
point, if the vg stats percentages differ from `mapping_metrics.tsv`, or if the BWA-MEM BAM does not have the primary records of the
giraffe BAMs. `chr20_table.py` makes the same checks on chromosomes 20-22. CMRG is not shown.

Supplementary File 3 holds the values of these figures and of Figure 7a; they are in the `summary.tsv` of each call set,
`benchmark_data.tsv`, `chr20_22_benchmark.tsv` (sheet Chromosomes 20-22) and the `graph_arms_*.tsv` tables (the script that wrote the Excel
file is not part of the release).

`chr20_22_strata_table.py` (not in `run_all.sh`) makes the tables of Figures 1 and 2 by region on chromosomes 20-22 in
`chr20_22_strata/` (`benchmark_data.tsv`, `graph_arms_data.tsv`, `graph_arms_roc.tsv`), from the `STRATA=1` runs of
`analyses/benchmark/bench_chr20.sh`; they are the sheet "Chromosomes 20-22 by region" of Supplementary File 3.

Colours: Paul Tol's muted scheme in Figure 1, S1 to S3; Okabe-Ito in Figure 2 and S5 to S7. Figure 1 and S3: Illumina dark colours
and circles, Aviti light colours and squares; in Figure 1 PASS calls solid, GATK before VQSR dashed; in S3 Illumina solid, Aviti
dashed. Figure 2 and S7: Illumina solid lines, Aviti dashed lines; every point a filled circle in the colour of its set-up; in 2e and
S7 a stem leads to each point. Figures 1 and 2 and S3 are drawn at the printed size, 183 mm wide (7 pt text, 8 pt titles).

Tools: python3 with matplotlib 3.6.3, numpy 1.26.4 and pillow, the versions of the run; other versions may change the bytes of the
images.
