# Figure 7: choice of the VQSR levels

`figure7_table.py` reads the hap.py summaries of the GATK call sets of HG002 and HG001 and writes `figure7_table.tsv`: recall,
precision, F1, false negatives and false positives of every call set, region, variant type and filter, with SNVs and indels
together as type ALL. It prints the comparisons quoted in the Methods, for example the level pair with the highest F1 in each
series. `figure7.py` draws Figure 7 from `figure7_table.tsv`.

## Inputs

The `summary.tsv` of each benchmark, copied into `data/` next to the scripts and named by the call set:

- `data/HG002/<Platform>_GATK_HaplotypeCaller_joint_vqsrbp_snp<S>_indel<I>.summary.tsv`: HG002 from `pipelines/gatk/vqsr`,
  benchmarked with `analyses/benchmark/happy.sh`.
- `data/HG001/HG001_<Platform>_GATK_HaplotypeCaller_<cohort|single>_vqsrbp_snp<S>_indel<I>.summary.tsv`, and for the cohort
  `..._vqsrbp_no_excesshet_snp99.9_indel99.9.summary.tsv`: hap.py against GIAB HG001 v4.2.1 (not in this release).

`<Platform>` is Illumina or Aviti, and `<S>/<I>` are the five level pairs 99.0/99.0, 99.7/99.0, 99.9/99.0, 99.9/99.7 and
99.9/99.9. Other sets in `data/` (the VQSR step of the joint-calling pipeline, `..._joint`, `..._single` and
`..._pipeline_vqsr`; DV and PPR on HG002) only enter the printed comparisons.

## Order

    python3 figure7_table.py > figure7_table.txt
    python3 figure7.py                    # figures/Figure7.pdf, .png and .tiff (600 dpi; DPI=300 for a smaller image)

`bash run_all.sh` runs both.

## Tools

Python 3 with matplotlib (3.11.2 for the published figure), numpy and pillow (the TIFF). Figure 7 is drawn at the printed size,
183 mm wide.
