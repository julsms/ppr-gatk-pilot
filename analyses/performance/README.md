# Computing performance

Execution time, CPU time and peak memory per sample of GATK and PPR (Supplementary Note 2, Supplementary Figure S28), from the
Nextflow execution traces of the per-sample runs. The structural-variant steps of PPR are left out; COMPRESS_INDEX_VCF (the
compression of the vg call VCF, about 13 s per sample) is counted. The QC steps of the GATK Aviti runs (MOSDEPTH, SAMTOOLS_STATS,
PICARD_ALIGNMENT_METRICS) are left out of CPU and memory; they ran in parallel and are part of the run duration. The GATK
per-sample runs also hold the per-sample GenomicsDB import (66 tasks; 204 = 6 + 3 x 66 tasks per run on Aviti, 336 = 6 + 5 x 66
on Illumina), and the time, CPU and memory include it.

## Inputs (all scripts run in this folder)

- `PPR/<platform>/<sample>_trace.txt`: the trace of the PPR run of the sample. Some of our runs held several samples, so the
  trace of each sample was rebuilt from the run traces and reports: one run per sample (the one with most of the 7 main steps),
  its COMPLETED and CACHED rows, the latest row per step, and a sample only with all 7 (KMC_COUNT, HAPLOTYPE_SAMPLING, VG_GIRAFFE,
  VG_VIEW, VG_INDEX_XG, VG_SURJECT, DEEPVARIANT). This kept 22 Illumina and 21 Aviti samples.
- `GATK/<platform>/individual/<sample>_trace.txt`: the trace of the per-sample GATK run. Leave out the traces of individual 23
  (`cpumin.sh` and `memoryRAM.sh` skip "23" only for GATK Illumina). `ppr_hosts.sh` pairs the GATK Aviti traces and times with
  the PPR samples by individual.
- The GATK time is the run duration in the Nextflow report (hours and minutes), typed into `average_time_persample.sh`.
- `ppr_v3_resources_per_sample.tsv`: per PPR sample, from the execution report of its run, the columns platform, sample, host,
  profile, cpus_heavy_step (the largest cpus of KMC_COUNT, HAPLOTYPE_SAMPLING, VG_GIRAFFE and VG_SURJECT), cpus_deepvariant and
  deepvariant_container.

## Order

    bash average_time_persample.sh        # total_time_results.csv: PPR, the sum of the task realtimes; GATK, the run duration
    bash cpumin.sh                        # cpu_minutes_results.csv: realtime x %cpu, summed over the tasks
    bash memoryRAM.sh                     # memory_usage_results.csv: the largest peak_rss of one task, and the sum
    bash ppr_hosts.sh                     # ppr_v3_matched_gatk_summary.tsv, ppr_v3_hosts.tsv, ppr_v3_hosts_summary.tsv

Each .csv table has one row per sample and an AVERAGE row per pipeline and platform. A GATK run is complete when its trace has the
most common number of tasks of the platform; `figures/performance/` leaves the other runs out of every metric. Column 24 of the
GATK Aviti traces is rss; `memoryRAM.sh` reads column 26 (peak_rss). The awk reads a realtime in ms (Nextflow writes it for tasks
under 1 s) as minutes, so check the realtime column first; our traces have none.

## Tools

bash, GNU awk, bc, column, Python 3 (standard library).
