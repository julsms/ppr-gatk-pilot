# Sequencing quality

Quality of the reads, the alignments and the calls of the cohort (Supplementary Table S2, Supplementary Note 1, Supplementary
File 2) and the kinship check of the 24 individuals, named by number (1-24). Individual 23 is left out of every summary (low yield on
Illumina), individual 6 only of the singleton summary.

## Inputs
Copy this folder to `<run>/qc` and run the scripts there, with the inputs below.
- `multiqc_illumina/<sample>/`: the Trim Galore reports and FastQC zips of the trimming before the pipelines (Trim Galore 0.6.10,
  cutadapt 4.4, `--fastqc`), and mosdepth and samtools flagstat of the GATK BAM from a QC step of the Illumina GATK runs whose
  command is not available; `bam_qc.sh` is a different call (the one used for the PPR BAMs).
- `multiqc_aviti/<sample>/`: the same reports, and the `qc/` folder of `pipelines/gatk/per_sample_aviti` (mosdepth, samtools stats).
- `bam_qc_ppr/<platform>_ppr/<sample>.*`: `bam_qc.sh` on the PPR BAM of each sample (projected onto GRCh38).
- In `PAPER_AGAIN_DIR` (default `..`, that is `<run>`): the outputs of `analyses/concordance` (`concordance/summary/all_counts.tsv`
  and the alleles of each individual in `concordance/prep/<set>_<platform>/`) and `uniquevariants/tables/r27_per_person.tsv`.
- For `kinship.sh`: per-sample VCFs (`vcfs/`) with the PASS non-reference genotypes of each sample. The GATK ones come from the annotated joint VCF (split
  VQSR 99.0) of an earlier joint calling of the 24 cohort genomes without HG002 (see `pipelines/gatk/README.md`), the PPR ones
  from an earlier PPR run.

## Order

    bash bam_qc.sh <PPR BAM> bam_qc_ppr/illumina_ppr/<sample>         # every sample, both platforms
    python3 generate_FASTqc_tables.py
    python3 generate_GATKbam_qc_tables.py
    python3 generate_PPRbam_qc_tables.py
    python3 generate_VCFqc_tables.py           # GATK, DV and PPR per individual, from the concordance alleles
    python3 coverage_gatk_ppr_correlation.py   # Pearson r of the GATK and PPR mean coverage (Supplementary Note 1)
    python3 generate_table1.py                 # Table1_manuscript.tsv
    bash kinship.sh GATK_illumina vcfs/GATK_illumina/*.GATK.filtered.vcf.gz   # also GATK_aviti, PPR_illumina, PPR_aviti

The generators write one per-individual table per platform and call set (`tables_*/`, the sheets of Supplementary File 2) and its
summary; the VCF tables come per individual from the concordance step. A release run leaves the VCF rows of
`Table1_manuscript.tsv` empty ("-"); `figures/qc/table1_revised.py` builds Supplementary Table S2 from its FASTQ and BAM rows only.

## Tools

mosdepth 0.3.6 and samtools 1.17 (biocontainers `mosdepth:0.3.6--hd299d5a_0`, `samtools:1.17--h00cdaf9_0`; samtools stats 1.19.2 in
`per_sample_aviti`), bcftools 1.20, PLINK 2.00a5.12 (`plink2:2.00a5.12--h4ac6f70_0`), Python 3.12.3 with pandas 2.1.4, numpy 1.26.4.
