# Graph decomposition: errors, effects and alignment metrics

The five set-ups of Figure 2 (BWA-MEM, GRCh38-only, HPRC default, HPRC d9, PPR) on HG002, one table set per platform.
`decompose.py` gives the total errors (FN + FP of SNVs and indels, PASS, as in Siren et al. 2024), the relative change of each
step and the share of each step. The changes against BWA-MEM in the abstract are the `rel_change_pct` of the steps
bwa_to_d9 (frequency-filtered graph) and total_vs_bwa (PPR) in `effects.tsv` (method `ga4gh_v5.0q`, region BB, type ALL). Supplementary Table S1 and the
changes of each set-up against BWA-MEM in the Results are made by `figures/benchmark/supp_table_graph_setups_errors.py` from the same
FN and FP (PASS, BB). `mapping_summary.py` gives the alignment metrics of the giraffe set-ups (Figure 2e and Supplementary Figure S7);
`bam_metrics.sh` on the BWA-MEM BAM gives those of BWA-MEM.

Run every step from one benchmark folder, the same as for `happy.sh` (or set `BENCHMARK_REVISION_DIR` to it). The scripts read:

| path in the benchmark folder | made by |
|---|---|
| `<call set>/summary.tsv`, `BENCH_DONE` | `<release>/analyses/benchmark/happy.sh` |
| `<arm>/mapping/vg_stats_a.txt`, `<arm>/work/<sample>_align.err` | `<release>/analyses/benchmark/graph_arms/giraffe_arm.sh` |
| `<Platform>_PPR_mapping_stats/mapping/` | `ppr_mapping_stats.sh`, from the GAM and giraffe log of the PPR pipeline |
| `comparisons/Mapping_metrics_graph_arms_<Platform>/bam/<arm>.tsv`, `bam_chr1_22XY/<arm>.tsv` (BWA-MEM) | `bam_metrics.sh` |
| `scripts/graph_arms/arms.tsv` | a copy of the file here |

## Order (in the benchmark folder; D=<release>/analyses/benchmark/decomposition)

```bash
mkdir -p scripts/graph_arms
cp $D/arms.tsv scripts/graph_arms/
bash $D/ppr_mapping_stats.sh <PPR GAM> <PPR giraffe log> Illumina_PPR_mapping_stats
bash $D/bam_metrics.sh <arm> <BAM> comparisons/Mapping_metrics_graph_arms_Illumina/bam   # three giraffe arms and PPR
python3 $D/decompose.py; python3 $D/mapping_summary.py                          # PLATFORM=aviti for Aviti
# BWA-MEM: the duplicate-marked BAM of the GATK pipeline (the DV input), on all contigs and on chr1-22, X, Y
bash $D/bam_metrics.sh Illumina_GATK_Deepvariant_markdup <BWA-MEM BAM> comparisons/Mapping_metrics_graph_arms_Illumina/bam
KEEP_CONTIGS=chr1_22XY bash $D/bam_metrics.sh Illumina_GATK_Deepvariant_markdup <BWA-MEM BAM> \
    comparisons/Mapping_metrics_graph_arms_Illumina/bam_chr1_22XY                 # the same for Aviti_..., in ..._Aviti/
```

Outputs: `comparisons/Decomposition_summary_with_Siren_<Platform>/` (`arm_metrics.tsv`, `effects.tsv`, `shares.tsv`,
`decomposition.md`; the folder name is that of the run) and `comparisons/Mapping_metrics_graph_arms_<Platform>/`
(`mapping_metrics.tsv`, `mapping_metrics.md`). BWA-MEM has no graph alignment, so it has only the BAM metrics. vg stats ran on the PPR GAM with 4 threads (Illumina) and
8 (Aviti); its output does not depend on the threads. `arms.tsv` names the machines where giraffe ran only by number.

Tools: python3, bcftools 1.20, tabix, samtools 1.19.2, vg 1.52.0.
