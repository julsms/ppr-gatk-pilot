# K-mer support of the variant calls (Merfin)

Merfin 1.1 checks each record of a call set against the k-mers of the reads of the same sample (`-filter` mode, defaults, no
`-peak` or `-prob`); a record has k-mer support when Merfin keeps it. The read database is one meryl count of the trimmed R1 and
R2 (k = 31) that keeps the k-mers seen more than once. The run counted the reads in chunks and summed them with `meryl union-sum`
to save disk; this gives the same counts. Direct mode (Merfin on the set itself) gives the rates of the paper.

## Inputs

- `sheet/samples.tsv` and `sheet/vcfs.tsv`: copy the templates into `<out>/sheet/` and list all 48 samples (23 individuals and
  HG002 on each platform), each with the eight sets of `config/merfin_sets.tsv`: the call sets prepared by `analyses/concordance`
  and the unique sets of `analyses/uniquevariants`; the HG002 sets come from `analyses/uniquevariants/hg002_precision.sh`.
- Give each individual its own reads (`analyses/genomescope/sample_pairs.csv` pairs the read sets of the two platforms): wrong
  reads change Figure 5, Supplementary Figures S19-S21 and Figure 6d.
- `GRCh38_chr.fa`: the chr-named GRCh38 of `pipelines/ppr` (the `REF` of `analyses/concordance`), with its `.fai` (md5
  `a368fcfa7065c744241ea5d2fa9ee431`).

## Order

Run from the folder that holds `concordance/`, `uniquevariants/` and `reads/`, with `kmers_validation/revision` as `<out>` (there
`figures/coverage/coverage_table.py` finds its table), and call the scripts by their path:

1. `bash ref_kmers.sh GRCh38_chr.fa <out>/ref/GRCh38.k31.meryl` (once).
2. `bash merfin_sample.sh <platform> <sample> <out> GRCh38_chr.fa <out>/ref/GRCh38.k31.meryl` per sample: `<set>.status.tsv.gz`
   (key, merfin 1/0) and `sample_counts.tsv` in `results/<platform>/<sample>/`; the read database in `db/` (20-30 GiB).
3. `bash tables.sh <out> uniquevariants`: `tables/` (per sample, cohort, paired tests) and `hg002/` (Merfin against GIAB). The
   cohort read databases can then be removed; `kstar_hg002.sh` needs the two HG002 ones.
4. `GENOMESCOPE=<path>/genomescope.R bash kstar_hg002.sh <platform> <out> GRCh38_chr.fa <out>/ref/GRCh38.k31.meryl uniquevariants`
   for HG002 on each platform: `hg002/kstar/<platform>/` (`kstar_2x2.tsv`, `kstar_summary.tsv`, `peak.tsv`, `verify.txt`).

The Python scripts take only the samples recorded in `SAMPLE_DONE` and `sample_state.tsv`, which `merfin_sample.sh` writes last.

## k* on HG002

The k-mers of the reference and allele sequence of each record (30 bases on each side) are counted in the reads and in GRCh38 with
meryl-lookup and turned into copy numbers as `merfin -dump` does, with the haploid peak of a GenomeScope 2.0 fit (`-k 31 -p 2`).
Supported calls are flagged for excess copies (k* >= 1) or repeated sequence; `verify.txt` gives the median copy number of the
supported GIAB true positives in single-copy context, and `merfin -dump` on about 2,000 records checks the counts.

## Tools

meryl 1.4.1 (r977) and meryl-lookup, Merfin 1.1 (r206), bcftools 1.20, htslib 1.21, GenomeScope 2.0 (`genomescope.R`, tag
v2.1.0), R 4.3.3, Python 3.12 with pandas 2.1.4, numpy 1.26.4, scipy 1.11.4. The Python code imports
`analyses/uniquevariants/lib/uvlib.py`, which reads `analyses/concordance/config/strata.tsv`. Threads (16 or 48) and meryl memory
(64 or 100 GB) do not change the counts.
