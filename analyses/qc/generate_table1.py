#!/usr/bin/env python3
"""The QC table of the cohort: the FASTQ, BAM and VCF summaries of the four generators in one table (Illumina and Aviti, GATK and
PPR), each cell "mean ± SD (min-max)" or "mean ± SD".
Input: the summary TSVs in tables_FATQC_illumina/, tables_FATSTQC_aviti/, tables_GATKbam_qc/, tables_PPRbam_qc/, tables_vcf_qc/
next to this script. Output: Table1_manuscript.tsv next to this script.
"""

import os
import csv
import sys


# ============================================================
# CONFIGURATION - File paths
# ============================================================
BASE_QC = "."
# Auto-detect: script directory is the QC base directory
BASE_QC = os.path.dirname(os.path.abspath(__file__))


SUMMARY_FILES = {
    # FASTQ level (platform-only, shared between pipelines)
    "fastq_illumina": os.path.join(BASE_QC, "tables_FATQC_illumina", "fastq_qc_summary_illumina.tsv"),
    "fastq_aviti": os.path.join(BASE_QC, "tables_FATSTQC_aviti", "fastq_qc_summary_aviti.tsv"),
    # BAM level (pipeline × platform)
    "bam_illumina_gatk": os.path.join(BASE_QC, "tables_GATKbam_qc", "bam_qc_summary_illumina_gatk.tsv"),
    "bam_aviti_gatk": os.path.join(BASE_QC, "tables_GATKbam_qc", "bam_qc_summary_aviti_gatk.tsv"),
    "bam_illumina_ppr": os.path.join(BASE_QC, "tables_PPRbam_qc", "bam_qc_summary_illumina_ppr.tsv"),
    "bam_aviti_ppr": os.path.join(BASE_QC, "tables_PPRbam_qc", "bam_qc_summary_aviti_ppr.tsv"),
    # VCF level (pipeline × platform)
    "vcf_gatk_illumina": os.path.join(BASE_QC, "tables_vcf_qc", "vcf_qc_summary_gatk_illumina.tsv"),
    "vcf_gatk_aviti": os.path.join(BASE_QC, "tables_vcf_qc", "vcf_qc_summary_gatk_aviti.tsv"),
    "vcf_ppr_illumina": os.path.join(BASE_QC, "tables_vcf_qc", "vcf_qc_summary_ppr_illumina.tsv"),
    "vcf_ppr_aviti": os.path.join(BASE_QC, "tables_vcf_qc", "vcf_qc_summary_ppr_aviti.tsv"),
}

OUTPUT_FILE = os.path.join(BASE_QC, "Table1_manuscript.tsv")


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def load_summary_tsv(filepath):
    """Load a summary TSV into a dict: metric_name -> {mean, sd, min, max, median, n}."""
    data = {}
    if not os.path.isfile(filepath):
        return data
    with open(filepath, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            metric = row.get("metric", "")
            try:
                data[metric] = {
                    "mean": float(row["mean"]) if row.get("mean", "NA") != "NA" else None,
                    "sd": float(row["sd"]) if row.get("sd", "NA") != "NA" else None,
                    "min": float(row["min"]) if row.get("min", "NA") != "NA" else None,
                    "max": float(row["max"]) if row.get("max", "NA") != "NA" else None,
                    "median": float(row["median"]) if row.get("median", "NA") != "NA" else None,
                    "n": int(row["n"]) if row.get("n", "0") not in ("NA", "0", "") else 0,
                }
            except (ValueError, KeyError):
                data[metric] = {"mean": None, "sd": None, "min": None, "max": None, "median": None, "n": 0}
    return data


def fmt_mean_sd(data_dict, metric, decimals=1, show_range=True):
    """Format a metric as 'mean ± SD (min-max)' or 'mean ± SD'."""
    if data_dict is None or metric not in data_dict:
        return "-"
    d = data_dict[metric]
    if d["mean"] is None:
        return "-"

    mean_str = f"{d['mean']:.{decimals}f}"
    sd_str = f"{d['sd']:.{decimals}f}" if d["sd"] is not None else "-"

    if show_range and d["min"] is not None and d["max"] is not None:
        min_str = f"{d['min']:.{decimals}f}"
        max_str = f"{d['max']:.{decimals}f}"
        return f"{mean_str} ± {sd_str} ({min_str}-{max_str})"
    else:
        return f"{mean_str} ± {sd_str}"


def fmt_mean_sd_int(data_dict, metric, show_range=True, divisor=1):
    """Format a metric with integer values (e.g., variant counts)."""
    if data_dict is None or metric not in data_dict:
        return "-"
    d = data_dict[metric]
    if d["mean"] is None:
        return "-"

    mean_val = d["mean"] / divisor
    sd_val = d["sd"] / divisor if d["sd"] is not None else 0

    if divisor > 1:
        mean_str = f"{mean_val:,.0f}"
        sd_str = f"{sd_val:,.0f}"
        if show_range and d["min"] is not None and d["max"] is not None:
            min_str = f"{d['min'] / divisor:,.0f}"
            max_str = f"{d['max'] / divisor:,.0f}"
            return f"{mean_str} ± {sd_str} ({min_str}-{max_str})"
        return f"{mean_str} ± {sd_str}"
    else:
        mean_str = f"{d['mean']:,.0f}"
        sd_str = f"{d['sd']:,.0f}" if d["sd"] is not None else "-"
        if show_range and d["min"] is not None and d["max"] is not None:
            min_str = f"{d['min']:,.0f}"
            max_str = f"{d['max']:,.0f}"
            return f"{mean_str} ± {sd_str} ({min_str}-{max_str})"
        return f"{mean_str} ± {sd_str}"


# ============================================================
# MAIN: Build Table 1
# ============================================================

def main():
    # Load all summary files
    summaries = {}
    missing_files = []
    for key, filepath in SUMMARY_FILES.items():
        summaries[key] = load_summary_tsv(filepath)
        if not summaries[key]:
            missing_files.append(key)

    if missing_files:
        print(f"WARNING: Could not load summary data for: {', '.join(missing_files)}")
        print("  These will appear as '-' in Table 1.")
        print(f"  (Expected paths under {BASE_QC})")
        print()

    # ---- Build Table 1 rows ----
    # FASTQ metrics: one column per platform (the same reads for both pipelines)
    # BAM and VCF metrics: Illumina GATK, Illumina PPR, Aviti GATK, Aviti PPR

    rows = []

    # ========================
    # SECTION: FASTQ-level QC
    # ========================
    rows.append(("SECTION", "FASTQ-level quality metrics", "", "", "", ""))

    # Header for FASTQ section (2 platforms)
    # Total reads (millions)
    rows.append((
        "Total reads (M)",
        fmt_mean_sd(summaries.get("fastq_illumina"), "total_reads_M", decimals=1),
        "",  # Illumina PPR = same as GATK (pre-alignment)
        fmt_mean_sd(summaries.get("fastq_aviti"), "total_reads_M", decimals=1),
        "",  # Aviti PPR = same as GATK
        ""
    ))

    # Adapter content (%)
    rows.append((
        "Reads with adapters (%)",
        fmt_mean_sd(summaries.get("fastq_illumina"), "adapters_pct", decimals=2),
        "",
        fmt_mean_sd(summaries.get("fastq_aviti"), "adapters_pct", decimals=2),
        "",
        ""
    ))

    # Quality-trimmed bases (%)
    rows.append((
        "Quality-trimmed bases (%)",
        fmt_mean_sd(summaries.get("fastq_illumina"), "quality_trimmed_pct", decimals=2),
        "",
        fmt_mean_sd(summaries.get("fastq_aviti"), "quality_trimmed_pct", decimals=2),
        "",
        ""
    ))

    # Total bases post-trimming (Gb)
    rows.append((
        "Total bases post-trim (Gb)",
        fmt_mean_sd(summaries.get("fastq_illumina"), "total_bases_Gb", decimals=1),
        "",
        fmt_mean_sd(summaries.get("fastq_aviti"), "total_bases_Gb", decimals=1),
        "",
        ""
    ))

    # Mean Phred quality
    rows.append((
        "Mean Phred quality score",
        fmt_mean_sd(summaries.get("fastq_illumina"), "mean_phred", decimals=1),
        "",
        fmt_mean_sd(summaries.get("fastq_aviti"), "mean_phred", decimals=1),
        "",
        ""
    ))

    # GC content (%)
    rows.append((
        "GC content (%)",
        fmt_mean_sd(summaries.get("fastq_illumina"), "pct_gc", decimals=1),
        "",
        fmt_mean_sd(summaries.get("fastq_aviti"), "pct_gc", decimals=1),
        "",
        ""
    ))

    # ========================
    # SECTION: BAM-level QC
    # ========================
    rows.append(("SECTION", "BAM-level alignment metrics", "", "", "", ""))

    # Mean coverage
    rows.append((
        "Mean coverage (×)",
        fmt_mean_sd(summaries.get("bam_illumina_gatk"), "mean_coverage", decimals=1),
        fmt_mean_sd(summaries.get("bam_illumina_ppr"), "mean_coverage", decimals=1),
        fmt_mean_sd(summaries.get("bam_aviti_gatk"), "mean_coverage", decimals=1),
        fmt_mean_sd(summaries.get("bam_aviti_ppr"), "mean_coverage", decimals=1),
        ""
    ))

    # % genome >= 10x
    rows.append((
        "Genome ≥10× (%)",
        fmt_mean_sd(summaries.get("bam_illumina_gatk"), "pct_10x", decimals=1),
        fmt_mean_sd(summaries.get("bam_illumina_ppr"), "pct_10x", decimals=1),
        fmt_mean_sd(summaries.get("bam_aviti_gatk"), "pct_10x", decimals=1),
        fmt_mean_sd(summaries.get("bam_aviti_ppr"), "pct_10x", decimals=1),
        ""
    ))

    # % genome >= 20x
    rows.append((
        "Genome ≥20× (%)",
        fmt_mean_sd(summaries.get("bam_illumina_gatk"), "pct_20x", decimals=1),
        fmt_mean_sd(summaries.get("bam_illumina_ppr"), "pct_20x", decimals=1),
        fmt_mean_sd(summaries.get("bam_aviti_gatk"), "pct_20x", decimals=1),
        fmt_mean_sd(summaries.get("bam_aviti_ppr"), "pct_20x", decimals=1),
        ""
    ))

    # % genome >= 30x
    rows.append((
        "Genome ≥30× (%)",
        fmt_mean_sd(summaries.get("bam_illumina_gatk"), "pct_30x", decimals=1),
        fmt_mean_sd(summaries.get("bam_illumina_ppr"), "pct_30x", decimals=1),
        fmt_mean_sd(summaries.get("bam_aviti_gatk"), "pct_30x", decimals=1),
        fmt_mean_sd(summaries.get("bam_aviti_ppr"), "pct_30x", decimals=1),
        ""
    ))

    # % mapped reads
    rows.append((
        "Mapped reads (%)",
        fmt_mean_sd(summaries.get("bam_illumina_gatk"), "pct_mapped", decimals=2),
        fmt_mean_sd(summaries.get("bam_illumina_ppr"), "pct_mapped", decimals=2),
        fmt_mean_sd(summaries.get("bam_aviti_gatk"), "pct_mapped", decimals=2),
        fmt_mean_sd(summaries.get("bam_aviti_ppr"), "pct_mapped", decimals=2),
        ""
    ))

    # % properly paired reads (samtools flagstat / samtools stats)
    rows.append((
        "Properly paired reads (%)",
        fmt_mean_sd(summaries.get("bam_illumina_gatk"), "pct_properly_paired", decimals=2),
        fmt_mean_sd(summaries.get("bam_illumina_ppr"), "pct_properly_paired", decimals=2),
        fmt_mean_sd(summaries.get("bam_aviti_gatk"), "pct_properly_paired", decimals=2),
        fmt_mean_sd(summaries.get("bam_aviti_ppr"), "pct_properly_paired", decimals=2),
        ""
    ))

    # ========================
    # SECTION: VCF-level QC
    # ========================
    rows.append(("SECTION", "VCF-level variant metrics (PASS-filtered)", "", "", "", ""))

    # Total variants
    rows.append((
        "Total variants",
        fmt_mean_sd_int(summaries.get("vcf_gatk_illumina"), "total_variants"),
        fmt_mean_sd_int(summaries.get("vcf_ppr_illumina"), "total_variants"),
        fmt_mean_sd_int(summaries.get("vcf_gatk_aviti"), "total_variants"),
        fmt_mean_sd_int(summaries.get("vcf_ppr_aviti"), "total_variants"),
        ""
    ))

    # SNV count
    rows.append((
        "SNV count",
        fmt_mean_sd_int(summaries.get("vcf_gatk_illumina"), "snv_count"),
        fmt_mean_sd_int(summaries.get("vcf_ppr_illumina"), "snv_count"),
        fmt_mean_sd_int(summaries.get("vcf_gatk_aviti"), "snv_count"),
        fmt_mean_sd_int(summaries.get("vcf_ppr_aviti"), "snv_count"),
        ""
    ))

    # INDEL count
    rows.append((
        "INDEL count",
        fmt_mean_sd_int(summaries.get("vcf_gatk_illumina"), "indel_count"),
        fmt_mean_sd_int(summaries.get("vcf_ppr_illumina"), "indel_count"),
        fmt_mean_sd_int(summaries.get("vcf_gatk_aviti"), "indel_count"),
        fmt_mean_sd_int(summaries.get("vcf_ppr_aviti"), "indel_count"),
        ""
    ))

    # Ti/Tv ratio
    rows.append((
        "Ti/Tv ratio",
        fmt_mean_sd(summaries.get("vcf_gatk_illumina"), "titv_ratio", decimals=2, show_range=False),
        fmt_mean_sd(summaries.get("vcf_ppr_illumina"), "titv_ratio", decimals=2, show_range=False),
        fmt_mean_sd(summaries.get("vcf_gatk_aviti"), "titv_ratio", decimals=2, show_range=False),
        fmt_mean_sd(summaries.get("vcf_ppr_aviti"), "titv_ratio", decimals=2, show_range=False),
        ""
    ))

    # Het/Hom ratio
    rows.append((
        "Het/Hom ratio",
        fmt_mean_sd(summaries.get("vcf_gatk_illumina"), "het_hom_ratio", decimals=2, show_range=False),
        fmt_mean_sd(summaries.get("vcf_ppr_illumina"), "het_hom_ratio", decimals=2, show_range=False),
        fmt_mean_sd(summaries.get("vcf_gatk_aviti"), "het_hom_ratio", decimals=2, show_range=False),
        fmt_mean_sd(summaries.get("vcf_ppr_aviti"), "het_hom_ratio", decimals=2, show_range=False),
        ""
    ))

    # Singleton count
    rows.append((
        "Singletons per sample†",
        fmt_mean_sd_int(summaries.get("vcf_gatk_illumina"), "singleton_count"),
        fmt_mean_sd_int(summaries.get("vcf_ppr_illumina"), "singleton_count"),
        fmt_mean_sd_int(summaries.get("vcf_gatk_aviti"), "singleton_count"),
        fmt_mean_sd_int(summaries.get("vcf_ppr_aviti"), "singleton_count"),
        ""
    ))

    # ========================
    # Write output TSV
    # ========================
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

    header = [
        "Metric",
        "Illumina GATK",
        "Illumina PPR",
        "Aviti GATK",
        "Aviti PPR",
        "Notes"
    ]

    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        for row in rows:
            if row[0] == "SECTION":
                # Write section header as a merged-like row
                writer.writerow([row[1], "", "", "", "", ""])
            else:
                writer.writerow(row)

    print(f"Table 1 written to: {OUTPUT_FILE}")
    print()

    # Also print to stdout for quick review
    print("=" * 120)
    print("TABLE 1: Sequencing quality assessment across platforms and variant calling methods")
    print("=" * 120)
    print(f"{'Metric':<28} {'Illumina GATK':<32} {'Illumina PPR':<32} {'Aviti GATK':<32} {'Aviti PPR':<32}")
    print("-" * 156)

    for row in rows:
        if row[0] == "SECTION":
            print(f"\n{'─' * 156}")
            print(f"  {row[1]}")
            print(f"{'─' * 156}")
        else:
            label = row[0]
            # For FASTQ rows, values are in col 1 (Illumina) and col 3 (Aviti)
            # PPR columns are empty for FASTQ
            col1 = row[1] if row[1] else "-"
            col2 = row[2] if row[2] else "-"
            col3 = row[3] if row[3] else "-"
            col4 = row[4] if row[4] else "-"
            print(f"{label:<28} {col1:<32} {col2:<32} {col3:<32} {col4:<32}")

    print()
    print("Notes:")
    print("  - Values are mean ± SD (min-max) unless otherwise noted.")
    print("  - FASTQ metrics are shared between GATK and PPR (pre-alignment processing is identical).")
    print("  - † Singleton counts exclude one of the two related samples (sample 6).")
    print("  - '-' marks a summary table that was not found.")
    print()


if __name__ == "__main__":
    main()
