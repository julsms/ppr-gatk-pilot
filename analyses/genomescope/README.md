# GenomeScope 2.0 on Illumina and Aviti reads

GenomeScope 2.0 fits a diploid model (k = 21) to the k-mer spectrum of the trimmed reads of each person on each platform. The
depth of a read set is the mosdepth mean coverage of its GATK alignments. The read set with the higher depth is first subsampled
with rasusa (seed 42) to the lower depth of the two platforms of the same person (`--target_bases` = floor(lower depth x 3.1e9)).
FastK counts the k-mers and `Histex -G` writes the histogram with a cumulative last bin. Then `statistical_analysis.py` removes
outliers (IQR x 3 on each of the four estimates, both platforms pooled), pairs the platforms by person and tests them with paired
Wilcoxon signed-rank tests.

## Order

1. For each person and platform, with the depth of the read set and the depth of the same person on the other platform (see
   below):
   `bash run_genomescope.sh <sample> <R1> <R2> <depth> <depth on the other platform> <gs>/genomescope_results_illumina`
   (`genomescope_results_aviti` for Aviti). Output: `<gs>/genomescope_results_<platform>/<sample>/`, with `<sample>_k21.histo`
   and `genomescope2_output/<sample>_summary.txt`.
2. `python3 statistical_analysis.py --illumina_dir <gs>/genomescope_results_illumina --aviti_dir <gs>/genomescope_results_aviti
   --sex_file <sex.tsv> --output_dir <gs>/statistical_results`. Output: `paired_data.csv` (21 persons in the paper),
   `paired_wilcoxon_results.csv`, `outliers_removed.csv`. `--sex_file` is required (two tab-separated columns without header:
   sample name, M or F); a path that does not exist skips the comparison between sexes, which the paper does not report.
3. Copy `sample_pairs.csv` into `<gs>`. The tests reported in the paper and the figures are made by `figures/genomescope/`
   (`GENOMESCOPE_DIR=<gs>`). The effect size of the paper is the matched-pairs rank-biserial r of `genomescope_stats.py`, not
   `effect_size_r` (z / sqrt(n)) of `statistical_analysis.py`.

Name each read set by its number: the paired table and `sample_pairs.csv` give the Illumina and Aviti read sets of the same
person. `statistical_analysis.py` also accepts the raw run names of the sequencing runs as folder and sample names.

The depths of step 1 are `mean_coverage` of `bam_qc_per_sample_<platform>_gatk.tsv` of `analyses/qc`; an Aviti read set takes
the row whose `source_label` is its number. The depth on the other platform is that of the same person (`sample_pairs.csv`). The
run read the depths from an earlier table with the same values.

## Tools

Nextflow 25.10.2 with Docker (`nextflow.config`: one machine, 60 cores and 100 GB; FastK with 35 threads and 70 GB). Containers
of `main.nf`: `community.wave.seqera.io/library/rasusa_fastk_seqtk:c8764ae5239ee72b` (rasusa, FastK v1.2 and Histex) and
`community.wave.seqera.io/library/genomescope2:2.1.0--c39bc08092274c2b` (GenomeScope 2.0). Python 3 with pandas, numpy and scipy.
