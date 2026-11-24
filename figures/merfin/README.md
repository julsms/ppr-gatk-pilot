# Figure 5 and the Merfin supplementary figures

K-mer support (Merfin) of the variant calls of GATK, BWA-MEM + DeepVariant (DV) and PPR: all calls of each pipeline, the variants
found by one pipeline and not by the other (GATK vs PPR, DV vs PPR) and the PPR-only variants against both; 23 individuals and
HG002 on Illumina and Aviti. The cohort rates are those of the direct mode (each set run through Merfin on its own).

## Input

The outputs of `analyses/kmers` in `KV_DIR` (default `kmers_validation/revision`): `tables/counts_long.tsv`,
`tables/completeness.tsv`, `hg002/merfin_giab_2x2.tsv`, `hg002/merfin_giab_summary.tsv` and
`hg002/kstar/<platform>/{kstar_2x2,kstar_summary}.tsv`, with their `TABLES_DONE`, `HG002_MERFIN_DONE` and `KSTAR_DONE` files.

## Order

`bash run_all.sh` runs, in this order:

| script | output |
|---|---|
| `merfin_table.py` | `mf_rates.tsv` (individuals numbered 1-23), `mf_hg002.tsv`, `mf_hg002_regions.tsv`, `mf_data.provenance.txt` |
| `figure5_combine.py` (with `figure5_panels.py`) | `figures/Figure5.{pdf,png,tiff}` |
| `supp_merfin_per_individual.py` | Supplementary Figure S19, k-mer support per individual |
| `supp_merfin_regions_cohort.py` | Supplementary Figure S20, k-mer support by region in the cohort |
| `supp_merfin_counts_by_region.py` | Supplementary Figure S21, unique variants by region, all and with k-mer support |
| `supp_merfin_hg002_giab_classes.py` | Supplementary Figure S23, k-mer support of the GIAB classes on HG002 |
| `supp_merfin_hg002_regions_aviti.py` | Supplementary Figure S22, Figure 5d for Aviti |

`mf_rates.tsv` numbers the 23 individuals 1-23 in numeric order (24 becomes 23). `merfin_table.py` stops when a check fails: 48
samples done, 23 individuals in every set, region and type, and the HG002 2x2 table equal to the two summaries. The figures go
to `figures/`: Figure 5 as PDF, PNG and TIFF (`DPI`, default 600), the supplementary figures as PDF and PNG at 600 dpi.

## Tools

Python 3 with numpy, matplotlib (the paper figures: matplotlib 3.6.3) and pillow (the TIFF of Figure 5, which is drawn at the
printed size, 183 mm wide).
