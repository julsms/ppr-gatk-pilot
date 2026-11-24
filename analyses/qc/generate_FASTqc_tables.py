#!/usr/bin/env python3
"""FASTQ QC per sample and per platform (mean, SD, min, max, median): reads, reads with adapters, quality-trimmed bases and bases
after trimming from the Trim Galore reports; mean Phred quality and GC from the FastQC zips (R1 and R2 averaged).
Input: one folder per sample in multiqc_illumina/ and multiqc_aviti/ (trimming_reports/ and the FastQC zips).
Output: tables_FATQC_illumina/ and tables_FATSTQC_aviti/ (fastq_qc_per_sample_<platform>.tsv, fastq_qc_summary_<platform>.tsv).
"""

import os
import re
import glob
import csv
import statistics
import zipfile


# ============================================================
# CONFIGURATION
# ============================================================
ILLUMINA_QC_DIR = "multiqc_illumina"
AVITI_QC_DIR = "multiqc_aviti"
OUTPUT_ILLUMINA = "tables_FATQC_illumina"
OUTPUT_AVITI = "tables_FATSTQC_aviti"

# Samples excluded from summary statistics (extremely low coverage)
# These are still included in per-sample tables but marked and excluded from summary
EXCLUDED_SAMPLES = {
    "illumina": ["23"],
    "aviti": ["23"],
}


# ============================================================
# PARSING FUNCTIONS
# ============================================================

def parse_trimming_report(filepath):
    """Parse TrimGalore trimming report .txt file.

    Extracts:
      - total_reads: number of reads processed
      - reads_with_adapters_pct: % reads containing adapter
      - quality_trimmed_pct: % bases removed by quality trimming
      - total_bp_written: total bases written after filtering
    """
    results = {
        "total_reads": None,
        "reads_with_adapters_pct": None,
        "quality_trimmed_pct": None,
        "total_bp_written": None,
    }
    if not os.path.isfile(filepath):
        return results

    with open(filepath, "r") as f:
        content = f.read()

    # Total reads processed
    m = re.search(r"Total reads processed:\s+([\d,]+)", content)
    if m:
        results["total_reads"] = int(m.group(1).replace(",", ""))

    # Reads with adapters
    m = re.search(r"Reads with adapters:\s+([\d,]+)\s+\(([\d.]+)%\)", content)
    if m:
        results["reads_with_adapters_pct"] = float(m.group(2))

    # Quality-trimmed bases
    m = re.search(r"Quality-trimmed:\s+([\d,]+) bp\s+\(([\d.]+)%\)", content)
    if m:
        results["quality_trimmed_pct"] = float(m.group(2))

    # Total written (filtered) - in bp
    m = re.search(r"Total written \(filtered\):\s+([\d,]+) bp", content)
    if m:
        results["total_bp_written"] = int(m.group(1).replace(",", ""))

    return results


def parse_fastqc_zip(filepath):
    """Parse FastQC zip to extract mean Phred quality score and %GC.

    Reads from fastqc_data.txt inside the zip:
      - 'Per sequence quality scores' section -> weighted mean Phred
      - 'Basic Statistics' section -> %GC
    """
    results = {"mean_phred": None, "pct_gc": None}
    if not os.path.isfile(filepath):
        return results

    try:
        with zipfile.ZipFile(filepath, "r") as zf:
            # Find fastqc_data.txt inside the zip
            data_file = None
            for name in zf.namelist():
                if name.endswith("fastqc_data.txt"):
                    data_file = name
                    break
            if data_file is None:
                return results

            with zf.open(data_file) as f:
                content = f.read().decode("utf-8")

        # Extract %GC from Basic Statistics section
        in_basic = False
        for line in content.split("\n"):
            if line.startswith(">>Basic Statistics"):
                in_basic = True
                continue
            if in_basic and line.startswith(">>END_MODULE"):
                break
            if in_basic and "%GC" in line:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    try:
                        results["pct_gc"] = float(parts[1])
                    except ValueError:
                        pass

        # Extract "Per sequence quality scores" section
        in_section = False
        total_count = 0
        weighted_sum = 0.0
        for line in content.split("\n"):
            if line.startswith(">>Per sequence quality scores"):
                in_section = True
                continue
            if in_section and line.startswith(">>END_MODULE"):
                break
            if in_section and not line.startswith("#") and line.strip():
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    try:
                        quality = float(parts[0])
                        count = float(parts[1])
                        weighted_sum += quality * count
                        total_count += count
                    except ValueError:
                        continue

        if total_count > 0:
            results["mean_phred"] = round(weighted_sum / total_count, 2)

    except (zipfile.BadZipFile, Exception):
        pass

    return results


# ============================================================
# FILE FINDING
# ============================================================

def find_trimming_reports(sample_dir):
    """Find R1 and R2 trimming reports in the sample directory."""
    paths = {"trim_r1": None, "trim_r2": None}

    # Check trimming_reports subdirectory first
    tr_dir = os.path.join(sample_dir, "trimming_reports")
    if not os.path.isdir(tr_dir):
        # Maybe they're directly in sample_dir
        tr_dir = sample_dir

    # Find R1
    r1_patterns = ["*r1*trimming_report.txt", "*R1*trimming_report.txt"]
    for pattern in r1_patterns:
        matches = glob.glob(os.path.join(tr_dir, pattern))
        if matches:
            paths["trim_r1"] = matches[0]
            break

    # Find R2
    r2_patterns = ["*r2*trimming_report.txt", "*R2*trimming_report.txt"]
    for pattern in r2_patterns:
        matches = glob.glob(os.path.join(tr_dir, pattern))
        if matches:
            paths["trim_r2"] = matches[0]
            break

    return paths


def find_fastqc_zips(sample_dir):
    """Find R1 and R2 FastQC zip files in the sample directory."""
    paths = {"fastqc_r1": None, "fastqc_r2": None}

    # Check fastqc_results/ subdirectory first, then fastqc/, then sample_dir
    fqc_dir = sample_dir
    for subdir_name in ["fastqc_results", "fastqc"]:
        candidate = os.path.join(sample_dir, subdir_name)
        if os.path.isdir(candidate):
            fqc_dir = candidate
            break

    # Find R1 zip - try multiple naming conventions
    r1_patterns = [
        "*_r1_*fastqc.zip", "*_R1_*fastqc.zip",
        "*r1*val_1*fastqc.zip", "*R1*val_1*fastqc.zip",
        "*val_1*fastqc.zip", "*_1_*fastqc.zip",
        "*r1*fastqc.zip", "*R1*fastqc.zip",
    ]
    for pattern in r1_patterns:
        matches = glob.glob(os.path.join(fqc_dir, pattern))
        if matches:
            paths["fastqc_r1"] = matches[0]
            break

    # Find R2 zip - try multiple naming conventions
    r2_patterns = [
        "*_r2_*fastqc.zip", "*_R2_*fastqc.zip",
        "*r2*val_2*fastqc.zip", "*R2*val_2*fastqc.zip",
        "*val_2*fastqc.zip", "*_2_*fastqc.zip",
        "*r2*fastqc.zip", "*R2*fastqc.zip",
    ]
    for pattern in r2_patterns:
        matches = glob.glob(os.path.join(fqc_dir, pattern))
        if matches:
            paths["fastqc_r2"] = matches[0]
            break

    return paths


# ============================================================
# SAMPLE QC EXTRACTION
# ============================================================

def extract_fastq_qc(sample_dir, sample_name):
    """Extract FASTQ-level QC metrics for a sample from trimming reports + FastQC."""
    trim_paths = find_trimming_reports(sample_dir)
    fqc_paths = find_fastqc_zips(sample_dir)
    row = {"sample": sample_name}

    trim_r1 = parse_trimming_report(trim_paths["trim_r1"]) if trim_paths["trim_r1"] else {}
    trim_r2 = parse_trimming_report(trim_paths["trim_r2"]) if trim_paths["trim_r2"] else {}

    # Total reads = R1 + R2 (in millions)
    r1_reads = trim_r1.get("total_reads")
    r2_reads = trim_r2.get("total_reads")
    if r1_reads is not None and r2_reads is not None:
        row["total_reads_M"] = round((r1_reads + r2_reads) / 1e6, 2)
    elif r1_reads is not None:
        row["total_reads_M"] = round(r1_reads / 1e6, 2)
    else:
        row["total_reads_M"] = None

    # Reads with adapters (%) = average of R1 and R2
    r1_adapt = trim_r1.get("reads_with_adapters_pct")
    r2_adapt = trim_r2.get("reads_with_adapters_pct")
    if r1_adapt is not None and r2_adapt is not None:
        row["adapters_pct"] = round((r1_adapt + r2_adapt) / 2, 2)
    elif r1_adapt is not None:
        row["adapters_pct"] = r1_adapt
    else:
        row["adapters_pct"] = None

    # Quality-trimmed bases (%) = average of R1 and R2
    r1_qt = trim_r1.get("quality_trimmed_pct")
    r2_qt = trim_r2.get("quality_trimmed_pct")
    if r1_qt is not None and r2_qt is not None:
        row["quality_trimmed_pct"] = round((r1_qt + r2_qt) / 2, 2)
    elif r1_qt is not None:
        row["quality_trimmed_pct"] = r1_qt
    else:
        row["quality_trimmed_pct"] = None

    # Total bases written (Gb) = R1 + R2
    r1_bp = trim_r1.get("total_bp_written")
    r2_bp = trim_r2.get("total_bp_written")
    if r1_bp is not None and r2_bp is not None:
        row["total_bases_Gb"] = round((r1_bp + r2_bp) / 1e9, 2)
    elif r1_bp is not None:
        row["total_bases_Gb"] = round(r1_bp / 1e9, 2)
    else:
        row["total_bases_Gb"] = None

    # Mean Phred quality (from FastQC) = average of R1 and R2
    fqc_r1 = parse_fastqc_zip(fqc_paths["fastqc_r1"]) if fqc_paths["fastqc_r1"] else {}
    fqc_r2 = parse_fastqc_zip(fqc_paths["fastqc_r2"]) if fqc_paths["fastqc_r2"] else {}

    r1_phred = fqc_r1.get("mean_phred")
    r2_phred = fqc_r2.get("mean_phred")
    if r1_phred is not None and r2_phred is not None:
        row["mean_phred"] = round((r1_phred + r2_phred) / 2, 2)
    elif r1_phred is not None:
        row["mean_phred"] = r1_phred
    else:
        row["mean_phred"] = None

    # %GC (from FastQC) = average of R1 and R2
    r1_gc = fqc_r1.get("pct_gc")
    r2_gc = fqc_r2.get("pct_gc")
    if r1_gc is not None and r2_gc is not None:
        row["pct_gc"] = round((r1_gc + r2_gc) / 2, 1)
    elif r1_gc is not None:
        row["pct_gc"] = r1_gc
    else:
        row["pct_gc"] = None

    return row


# ============================================================
# STATISTICS
# ============================================================

def compute_summary_stats(rows, columns):
    """Compute mean, SD, min, max, median for numeric columns."""
    summary = {}
    for col in columns:
        values = [r[col] for r in rows if r[col] is not None]
        if not values:
            summary[col] = {"mean": None, "sd": None, "min": None, "max": None, "median": None, "n": 0}
            continue
        summary[col] = {
            "mean": round(statistics.mean(values), 2),
            "sd": round(statistics.stdev(values), 2) if len(values) > 1 else 0,
            "min": round(min(values), 2),
            "max": round(max(values), 2),
            "median": round(statistics.median(values), 2),
            "n": len(values),
        }
    return summary


# ============================================================
# MAIN PROCESSING
# ============================================================

def process_platform(qc_dir, output_dir, platform):
    """Process all samples for a platform, write per-sample TSV + summary."""
    os.makedirs(output_dir, exist_ok=True)

    if not os.path.isdir(qc_dir):
        print(f"  ERROR: QC directory not found: {qc_dir}")
        return

    sample_dirs = sorted([
        d for d in os.listdir(qc_dir)
        if os.path.isdir(os.path.join(qc_dir, d))
    ])

    if not sample_dirs:
        print(f"  WARNING: No sample directories found in {qc_dir}")
        return

    # Output columns (FASTQ-level)
    columns = [
        "total_reads_M",
        "adapters_pct",
        "quality_trimmed_pct",
        "total_bases_Gb",
        "mean_phred",
        "pct_gc",
    ]

    # Determine which samples to exclude from summary
    excluded = set(EXCLUDED_SAMPLES.get(platform, []))

    all_rows = []
    for sample_name in sample_dirs:
        sample_dir = os.path.join(qc_dir, sample_name)
        row = extract_fastq_qc(sample_dir, sample_name)
        row["excluded"] = sample_name in excluded
        all_rows.append(row)

    # Rows included in summary (excluding low-coverage samples)
    included_rows = [r for r in all_rows if not r["excluded"]]

    if excluded:
        print(f"  Excluded from summary: {sorted(excluded)} (low coverage)")
    print(f"  Included in summary: {len(included_rows)} / {len(all_rows)} samples")

    # Write per-sample table (ALL samples, with exclusion note column)
    per_sample_file = os.path.join(output_dir, f"fastq_qc_per_sample_{platform}.tsv")
    with open(per_sample_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample"] + columns + ["excluded"], delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in all_rows:
            out_row = {k: (v if v is not None else "NA") for k, v in row.items()}
            out_row["excluded"] = "yes" if row.get("excluded") else "no"
            writer.writerow(out_row)

    print(f"  Per-sample table: {per_sample_file}")
    print(f"  Total samples: {len(all_rows)}")

    # Compute and write summary (only included samples)
    summary = compute_summary_stats(included_rows, columns)
    summary_file = os.path.join(output_dir, f"fastq_qc_summary_{platform}.tsv")
    with open(summary_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["metric", "mean", "sd", "min", "max", "median", "n"])
        for col in columns:
            s = summary[col]
            writer.writerow([
                col,
                s["mean"] if s["mean"] is not None else "NA",
                s["sd"] if s["sd"] is not None else "NA",
                s["min"] if s["min"] is not None else "NA",
                s["max"] if s["max"] is not None else "NA",
                s["median"] if s["median"] is not None else "NA",
                s["n"],
            ])

    print(f"  Summary table: {summary_file}")

    # Print summary to stdout
    print(f"\n  === {platform.upper()} FASTQ QC Summary ===")
    print(f"  {'Metric':<22} {'Mean':>10} {'SD':>8} {'Min':>10} {'Max':>10} {'Median':>10} {'N':>4}")
    print(f"  {'-'*22} {'-'*10} {'-'*8} {'-'*10} {'-'*10} {'-'*10} {'-'*4}")
    for col in columns:
        s = summary[col]
        vals = [
            f"{s['mean']}" if s['mean'] is not None else "NA",
            f"{s['sd']}" if s['sd'] is not None else "NA",
            f"{s['min']}" if s['min'] is not None else "NA",
            f"{s['max']}" if s['max'] is not None else "NA",
            f"{s['median']}" if s['median'] is not None else "NA",
        ]
        print(f"  {col:<22} {vals[0]:>10} {vals[1]:>8} {vals[2]:>10} {vals[3]:>10} {vals[4]:>10} {s['n']:>4}")


def main():
    print("=" * 60)
    print("FASTQ QC Table Generation (Trimming Reports + FastQC)")
    print("=" * 60)

    print(f"\nProcessing ILLUMINA...")
    process_platform(ILLUMINA_QC_DIR, OUTPUT_ILLUMINA, "illumina")

    print(f"\nProcessing AVITI...")
    process_platform(AVITI_QC_DIR, OUTPUT_AVITI, "aviti")

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)


if __name__ == "__main__":
    main()
