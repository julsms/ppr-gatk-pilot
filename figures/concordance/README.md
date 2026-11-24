# Concordance figures

Figure 3, Supplementary Figures S8-S12 and the genome and BB rows of Supplementary File 4, made from `summary/all_counts.tsv` of
`analyses/concordance`: 23 individuals, the nine call-set pairs, the genome and the 10 regions, all variants, SNVs and indels. The
rows of the other nine regions of Supplementary File 4 and its exact-key columns (L0_*) are rows and columns of `all_counts.tsv`;
the workbook itself is built outside this release.

## Order

`CONCORDANCE_DIR=<absolute path of the concordance output folder> bash run_all.sh` runs the scripts in this order.

| Script | Output |
|---|---|
| `concordance_table.py` | `concordance_data.tsv` (concordance per individual) and `concordance_tests.tsv` |
| `figure3_panels.py`, `figure3_combine.py` | `figures/Figure3a_pairs`, `Figure3b_same_pipeline_vs_same_platform`, `Figure3c_across_platforms`, `Figure3` |
| `supp_concordance_crossplatform_tests.py` | Supplementary Figure S12 |
| `supp_concordance_caller_vs_reference.py` | Supplementary Figure S10 |
| `supp_concordance_all_pairs.py` | Supplementary Figure S9 |
| `supp_concordance_across_platforms_lines.py` | Supplementary Figure S11 |
| `supp_concordance_heatmap.py` | Supplementary Figure S8 |
| `supp_table_concordance_counts.py` | `SuppTable_concordance_counts.tsv` (the genome and BB rows of Supplementary File 4) |

Tables go into this folder, figures into `figures/` (PDF and PNG at 600 dpi; Figure 3 also TIFF). `concordance_data.tsv` numbers the
23 individuals 1-23 in the order of their `all_counts.tsv` numbers sorted as text (1, 10, ..., 19, 2, 20, 21, 22, 24, 3, ..., 9);
`SuppTable_concordance_counts.tsv` keeps those numbers (1-24 without 23), as Supplementary File 4 does.

## Tests

Paired Wilcoxon signed-rank test per region, n = 23: scipy `wilcoxon` on the non-zero differences, with continuity correction and
the normal approximation. Holm correction over the 10 regions for each comparison, variant type and metric; the genome row keeps its
raw p. Effect sizes: r = |z| / sqrt(n) and the matched-pairs rank-biserial correlation.

`posthoc_long.tsv` holds the tests of the original run, made with Holm over 12 regions by `07_stats_figures.py` (not in this
release): one row per comparison, region, type and metric (J_gt, J_al, and L0_gt_conc, the exact-key genotype concordance), with n,
the means and their difference in percentage points, V, raw p, both effect sizes, the Holm family size and adjusted p, the role of
the row (stratum, genome summary or nested in genome) and in_submitted, a True/False column that no release script reads.
`concordance_table.py` checks its means, n, V, raw p and effect sizes against this table and stops at a difference;
`run_all.sh` copies it into `$CONCORDANCE_DIR/stats/` first.

## Tools

Python 3.12.3 with numpy 1.26.4, scipy 1.11.4 and matplotlib 3.6.3; pillow for the TIFF of Figure 3, which is drawn at the printed
size, 183 mm wide.
