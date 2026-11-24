# GenomeScope supplementary figures

GenomeScope 2.0 estimates from the depth-matched reads of each individual on Illumina and Aviti (`analyses/genomescope`).
`genomescope_stats.py` tests the platforms with paired Wilcoxon signed-rank tests on the paired table of that analysis (21
individuals after its outlier rule), Holm-corrected over the four estimates, with the matched-pairs rank-biserial correlation as
effect size, and gives the Spearman correlations of the heterozygosity with the other estimates. It also checks that W and p equal
those of `paired_wilcoxon_results.csv`. The effect size of the paper is this rank-biserial r; the column `statistical_analysis_r`
only repeats `effect_size_r` (z / sqrt(n)) of `statistical_analysis.py`. The k-mer spectra are drawn for all 23 individuals.

## Inputs

- `GENOMESCOPE_DIR` (default `genomescope`): `statistical_results/` (`paired_data.csv`, `paired_wilcoxon_results.csv`,
  `outliers_removed.csv`), `sample_pairs.csv`, and `genomescope_results_{illumina,aviti}/<sample>/<sample>_k21.histo`;
- `PAPER_AGAIN_DIR` (default: the current folder): `uniquevariants/tables/r27_per_person.tsv` of the unique-variant analysis, for
  the heterozygous alleles in the call sets (`het_from_calls.py`, `het_platform_test.py`).

## Order

`bash run_all.sh` runs:

| script | output |
|---|---|
| `genomescope_stats.py` | `genomescope_stats.tsv`, `genomescope_correlations.tsv`, `genomescope_per_individual.tsv` (1-21) |
| `het_from_calls.py` | `het_from_calls.tsv`: heterozygous alleles per individual in the call sets, for the Results text |
| `het_platform_test.py` | `het_platform_test.tsv`: heterozygous rate and het/hom of the call sets, Illumina vs Aviti (paired Wilcoxon, Holm over the three pipelines; the Results give the 21 individuals of the GenomeScope comparison) |
| `supp_kmer_spectra.py` | Supplementary Figure S25, k-mer spectra of the paired samples |
| `supp_genomescope_estimates.py` (with `genomescope_panels.py`) | Supplementary Figure S26, the four estimates |
| `supp_genomescope_het_vs_length.py` (with `genomescope_panels.py`) | Supplementary Figure S27, heterozygosity against genome length |

Tables are written next to the scripts and figures to `figures/` (PNG at 600 dpi, environment variable `DPI`).
`genomescope_per_individual.tsv` numbers the 21 individuals kept by the outlier rule 1-21 in the order of `paired_data.csv`, which
is the text order of the Illumina folder names. `genomescope_stats.py` also accepts the raw run names of the sequencing runs as
sample names.

## Tools

Python 3 with numpy, scipy and matplotlib (the paper figures: matplotlib 3.6.3 for the spectra, 3.11.2 for S26 and S27).
