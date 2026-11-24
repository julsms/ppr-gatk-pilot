# Supplementary Figure S28: computing performance

`performance_table.py` reads the tables of `analyses/performance` and writes `performance_data.tsv` (one value per sample and
metric, the samples numbered) and `performance_summary.tsv` (n, mean, SD and range per pipeline, platform and metric, and the PPR
and GATK execution times on the same individuals). These are the values of Supplementary Note 2. GATK runs with an incomplete
trace (`gatk_run_completeness.tsv`, not made by this release: one row per run with platform, sample and complete yes/no) are
left out of every metric. It stops with an error if the mean of a table differs from its
AVERAGE row. `supp_performance.py` draws the four panels from `performance_data.tsv`: a execution time, b CPU time, c peak memory
of the most memory-intensive process, d peak memory summed over processes (bars: mean; whiskers: SD).

## Inputs

From `analyses/performance`: `total_time_results.csv`, `cpu_minutes_results.csv`, `memory_usage_results.csv`,
`gatk_run_completeness.tsv` and `ppr_v3_matched_gatk_summary.tsv`.

## Order

    TIMELINE_DIR=../../analyses/performance python3 performance_table.py
    python3 supp_performance.py           # figures/SuppFig_performance.png (600 dpi) and .pdf

`TIMELINE_DIR=<folder> bash run_all.sh` runs both. They write next to themselves; `performance_data.provenance.txt` lists the md5
of the inputs and outputs.

## Tools

Python 3 with matplotlib and numpy.
