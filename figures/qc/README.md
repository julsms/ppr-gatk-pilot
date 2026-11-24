# Supplementary Table S2: sequencing, alignment and variant-call quality

`table1_revised.py` builds Supplementary Table S2 with six columns, GATK, DV and PPR on Illumina and on Aviti. Each cell is the
mean ± SD over the 23 individuals with the range in parentheses (only mean ± SD for Ti/Tv and het/hom).

- FASTQ and BAM rows: `Table1_manuscript.tsv` of `analyses/qc`. The FASTQ rows are given once per platform, because the three
  pipelines used the same reads. DV used the GATK alignments, so it shares the GATK BAM rows.
- VCF rows: the PASS variant alleles carried by each individual after splitting and normalization, as in the concordance
  analysis: the totals of `concordance/summary/all_counts.tsv`, and Ti/Tv and het/hom from the four concordance classes of
  `uniquevariants/tables/r27_per_person.tsv`. Counts are in millions.

The script stops with an error if a call set does not have 23 individuals, if the total of a call set differs between the pairs it
is in, or if the four classes do not add up to the total of an individual.

## Order

    PAPER_AGAIN_DIR=<folder with qc/, concordance/summary/ and uniquevariants/tables/> python3 table1_revised.py

It writes `Table1_revised.tsv` (20 rows) and `Table1_revised.provenance.txt` (md5 of the inputs and of the table) next to itself.

## Tools

Python 3 (standard library).
