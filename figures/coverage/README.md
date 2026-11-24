# Figure 6: sequencing depth, concordance and k-mer support

Depth is the mosdepth mean coverage of the GATK alignments (BWA-MEM, duplicates marked) of each individual. `coverage_table.py`
correlates it with the genotype concordance of the call-set pairs and with the k-mer support of the call sets, for the 23
individuals on each platform: Spearman rho (and Pearson r), Holm-corrected within each panel family of Figure 6 (genome, all
variants) and over the 10 regions (regions). The genome rows for SNVs and indels and the correlation of the Illumina and Aviti
depths keep their raw p; `holm_by_variant_type.py` corrects the SNV and indel rows as the Results report them. Between platforms
the depth is the lower of the two depths, and also their difference. The figures add Theil-Sen guide lines.

## Inputs

Under `PAPER_AGAIN_DIR` (default: the current folder):

- `qc/tables_GATKbam_qc/bam_qc_per_sample_{illumina,aviti}_gatk.tsv`: depth per sample (`analyses/qc`). The Aviti table must be
  numbered by individual, with the label of the source file in a `source_label` column;
- `concordance/summary/all_counts.tsv`: genotype concordance per individual, pair, region and type (`analyses/concordance`);
- `kmers_validation/revision/tables/per_sample_genome.tsv`: k-mer support per sample and call set (`analyses/kmers`).

`coverage_table.py` also accepts the raw run names of the sequencing runs as QC sample names.

## Order

`bash run_all.sh` runs:

1. `coverage_table.py`: `coverage_data.tsv`, `coverage_correlations.tsv` and `coverage_table.provenance.txt`, next to the
   script. `coverage_data.tsv` numbers the 23 individuals 1-23 in the order of their numbers sorted as text (1, 10, 11, ..., 19,
   2, 20, 21, 22, 24, 3, ..., 9), as `figures/concordance` does;
2. `figure6.py` (with `figure6_panels.py`): `figures/Figure6.{pdf,png,tiff}`;
3. `supp_depth_by_region.py`: `figures/SuppFig_depth_by_region.{pdf,png}` (Supplementary Figure S24);
4. `holm_by_variant_type.py`: `coverage_correlations_by_type.tsv`, the genome rows of SNVs and indels with Holm over the three
   pipelines (depth difference) or over the three pipeline pairs of a platform (depth).

The correlation code sits here, next to the figure scripts, because they read its tables from this folder.

## Tools

Python 3 with numpy, scipy, matplotlib (the paper figures: matplotlib 3.6.3) and pillow (the TIFF of Figure 6, which is drawn at the
printed size, 183 mm wide).
