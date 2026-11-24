#!/usr/bin/env python3
"""BAM QC of the GATK alignments per sample and per platform (mean, SD, min, max, median): mean coverage and % of the genome at
10x, 20x and 30x (mosdepth summary and global.dist), % mapped and % properly paired reads (samtools flagstat, or samtools stats).
Input: one folder per sample in multiqc_illumina/ (flat files) and multiqc_aviti/ (mosdepth/ and samtools/ subfolders).
Output: tables_GATKbam_qc/bam_qc_per_sample_<platform>_gatk.tsv and bam_qc_summary_<platform>_gatk.tsv.
"""

import os
import re
import glob
import csv
import statistics

# ============================================================
# CONFIGURATION
# ============================================================
ILLUMINA_GATK_DIR = "multiqc_illumina"
AVITI_GATK_DIR = "multiqc_aviti"
OUTPUT_DIR = "tables_GATKbam_qc"

# Samples excluded from summary (extremely low coverage)
EXCLUDED_SAMPLES = {
    "illumina_gatk": ["23"],
    "aviti_gatk": ["23"],
}

# ============================================================
# PARSING FUNCTIONS
# ============================================================

def parse_mosdepth_summary(filepath):
    """Parse mosdepth summary file to get mean coverage.
    
    Format: chrom\tlength\tbases\tmean\tmin\tmax
    We want the 'total' or 'total_region' row.
    """
    if not os.path.isfile(filepath):
        return None
    
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                # Look for total or total_region row
                if parts[0] in ("total", "total_region"):
                    try:
                        return float(parts[3])
                    except (ValueError, IndexError):
                        pass
    return None


def parse_mosdepth_global_dist(filepath):
    """Parse mosdepth global.dist.txt to get % genome at thresholds.
    
    Format: chrom\tcoverage\tfraction_at_or_above
    We look at 'total' rows for coverage = 10, 20, 30.
    The fraction is cumulative (proportion of bases >= that coverage).
    """
    results = {"pct_10x": None, "pct_20x": None, "pct_30x": None}
    if not os.path.isfile(filepath):
        return results
    
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3 and parts[0] == "total":
                try:
                    cov = int(parts[1])
                    frac = float(parts[2])
                    if cov == 10:
                        results["pct_10x"] = round(frac * 100, 2)
                    elif cov == 20:
                        results["pct_20x"] = round(frac * 100, 2)
                    elif cov == 30:
                        results["pct_30x"] = round(frac * 100, 2)
                except (ValueError, IndexError):
                    pass
    return results


def parse_mosdepth_thresholds(filepath):
    """Parse mosdepth thresholds file (if --thresholds was used).
    
    Format: #chrom\tstart\tend\t10X\t20X\t30X
    Or summary: chrom\t10X\t20X\t30X (proportions)
    
    This is an alternative to global.dist for getting coverage thresholds.
    """
    # Typically the thresholds bed.gz is harder to parse; 
    # we'll rely on global.dist instead.
    return {"pct_10x": None, "pct_20x": None, "pct_30x": None}


def parse_flagstat(filepath):
    """Parse samtools flagstat output to get % mapped reads.
    
    Example line: "123456789 + 0 mapped (99.45% : N/A)"
    """
    if not os.path.isfile(filepath):
        return None
    
    with open(filepath, "r") as f:
        content = f.read()
    
    # Look for mapped line with percentage
    m = re.search(r"(\d+)\s+\+\s+\d+\s+mapped\s+\((\d+\.?\d*)%", content)
    if m:
        return float(m.group(2))
    return None


def parse_samtools_stats_mapped(filepath):
    """Parse samtools stats output to get % mapped reads.
    
    Look for 'SN\treads mapped:' and 'SN\traw total sequences:' lines.
    """
    if not os.path.isfile(filepath):
        return None
    
    total_reads = None
    mapped_reads = None
    
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("SN\traw total sequences:"):
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    try:
                        total_reads = int(parts[2])
                    except ValueError:
                        pass
            elif line.startswith("SN\treads mapped:"):
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    try:
                        mapped_reads = int(parts[2])
                    except ValueError:
                        pass
    
    if total_reads and mapped_reads:
        return round((mapped_reads / total_reads) * 100, 2)
    return None



def parse_flagstat_properly_paired(filepath):
    """Get % properly paired from flagstat ("... properly paired (93.75% : N/A)")."""
    if not os.path.isfile(filepath):
        return None
    with open(filepath, "r") as f:
        content = f.read()
    m = re.search(r"(\d+)\s+\+\s+\d+\s+properly paired\s+\((\d+\.?\d*)%", content)
    if m:
        return float(m.group(2))
    return None


def parse_samtools_stats_properly_paired(filepath):
    """% properly paired from samtools stats: 'SN reads properly paired' / 'SN raw total sequences' * 100
    (all reads are paired-end, so this equals flagstat's properly-paired percentage)."""
    if not os.path.isfile(filepath):
        return None
    total_reads = None
    proper = None
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("SN\traw total sequences:"):
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    try:
                        total_reads = int(parts[2])
                    except ValueError:
                        pass
            elif line.startswith("SN\treads properly paired:"):
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    try:
                        proper = int(parts[2])
                    except ValueError:
                        pass
    if total_reads and proper is not None:
        return round((proper / total_reads) * 100, 2)
    return None

# ============================================================
# FILE FINDING (supports both layouts)
# ============================================================

def find_mosdepth_files(sample_dir, sample_name):
    """Find mosdepth summary and global dist files."""
    paths = {"summary": None, "global_dist": None}
    
    # Try Layout B first (subdirectory)
    mosdepth_dir = os.path.join(sample_dir, "mosdepth")
    if os.path.isdir(mosdepth_dir):
        search_dir = mosdepth_dir
    else:
        search_dir = sample_dir
    
    # Find summary
    patterns_summary = [
        f"{sample_name}.mosdepth.summary.txt",
        f"{sample_name}*.mosdepth.summary.txt",
    ]
    for pattern in patterns_summary:
        matches = glob.glob(os.path.join(search_dir, pattern))
        if matches:
            paths["summary"] = matches[0]
            break
    
    # If not found in subdir, try flat
    if paths["summary"] is None:
        for pattern in patterns_summary:
            matches = glob.glob(os.path.join(sample_dir, pattern))
            if matches:
                paths["summary"] = matches[0]
                break
    
    # Find global.dist
    patterns_dist = [
        f"{sample_name}.mosdepth.global.dist.txt",
        f"{sample_name}*.mosdepth.global.dist.txt",
    ]
    for pattern in patterns_dist:
        matches = glob.glob(os.path.join(search_dir, pattern))
        if matches:
            paths["global_dist"] = matches[0]
            break
    
    if paths["global_dist"] is None:
        for pattern in patterns_dist:
            matches = glob.glob(os.path.join(sample_dir, pattern))
            if matches:
                paths["global_dist"] = matches[0]
                break
    
    return paths


def find_mapped_reads_file(sample_dir, sample_name):
    """Find flagstat or samtools stats file for % mapped."""
    
    # Try flagstat first (Layout A - flat)
    patterns_flagstat = [
        f"{sample_name}_flagstat.txt",
        f"{sample_name}*flagstat*.txt",
    ]
    for pattern in patterns_flagstat:
        matches = glob.glob(os.path.join(sample_dir, pattern))
        if matches:
            return ("flagstat", matches[0])
    
    # Try samtools subdir (Layout B)
    samtools_dir = os.path.join(sample_dir, "samtools")
    if os.path.isdir(samtools_dir):
        patterns_stats = [
            f"{sample_name}_samtools_stats.txt",
            f"{sample_name}*stats*.txt",
        ]
        for pattern in patterns_stats:
            matches = glob.glob(os.path.join(samtools_dir, pattern))
            if matches:
                return ("stats", matches[0])
    
    # Try stats in flat dir
    patterns_stats_flat = [
        f"{sample_name}_stats.txt",
        f"{sample_name}*_stats.txt",
    ]
    for pattern in patterns_stats_flat:
        matches = glob.glob(os.path.join(sample_dir, pattern))
        if matches:
            return ("stats", matches[0])
    
    return (None, None)


# ============================================================
# SAMPLE QC EXTRACTION
# ============================================================

def extract_bam_qc(sample_dir, sample_name):
    """Extract BAM-level QC metrics for a sample."""
    row = {"sample": sample_name}
    
    # --- Coverage from mosdepth ---
    mosdepth_paths = find_mosdepth_files(sample_dir, sample_name)
    
    # Mean coverage
    if mosdepth_paths["summary"]:
        row["mean_coverage"] = parse_mosdepth_summary(mosdepth_paths["summary"])
    else:
        row["mean_coverage"] = None
    
    # % genome at thresholds
    if mosdepth_paths["global_dist"]:
        thresholds = parse_mosdepth_global_dist(mosdepth_paths["global_dist"])
        row["pct_10x"] = thresholds["pct_10x"]
        row["pct_20x"] = thresholds["pct_20x"]
        row["pct_30x"] = thresholds["pct_30x"]
    else:
        row["pct_10x"] = None
        row["pct_20x"] = None
        row["pct_30x"] = None
    
    # --- % mapped reads ---
    file_type, filepath = find_mapped_reads_file(sample_dir, sample_name)
    if file_type == "flagstat":
        row["pct_mapped"] = parse_flagstat(filepath)
    elif file_type == "stats":
        row["pct_mapped"] = parse_samtools_stats_mapped(filepath)
    else:
        row["pct_mapped"] = None

    # --- % properly paired reads ---
    if file_type == "flagstat":
        row["pct_properly_paired"] = parse_flagstat_properly_paired(filepath)
    elif file_type == "stats":
        row["pct_properly_paired"] = parse_samtools_stats_properly_paired(filepath)
    else:
        row["pct_properly_paired"] = None
    
    return row


# ============================================================
# STATISTICS
# ============================================================

def compute_summary_stats(rows, columns):
    """Compute mean, SD, min, max, median for numeric columns."""
    summary = {}
    for col in columns:
        values = [r[col] for r in rows if r.get(col) is not None]
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

def process_platform(qc_dir, output_prefix, platform_key):
    """Process all samples for a platform/pipeline, write per-sample TSV + summary."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
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
    
    columns = ["mean_coverage", "pct_10x", "pct_20x", "pct_30x", "pct_mapped", "pct_properly_paired"]
    excluded = set(EXCLUDED_SAMPLES.get(platform_key, []))
    
    all_rows = []
    missing_data = []
    
    for sample_name in sample_dirs:
        sample_dir = os.path.join(qc_dir, sample_name)
        row = extract_bam_qc(sample_dir, sample_name)
        row["excluded"] = sample_name in excluded
        all_rows.append(row)
        
        # Track missing data
        missing = [col for col in columns if row.get(col) is None]
        if missing:
            missing_data.append((sample_name, missing))
    
    # Report missing data
    if missing_data:
        print(f"  WARNING: Missing data for {len(missing_data)} sample(s):")
        for sname, mcols in missing_data[:5]:
            print(f"    {sname}: missing {', '.join(mcols)}")
        if len(missing_data) > 5:
            print(f"    ... and {len(missing_data) - 5} more")
    
    # Rows included in summary
    included_rows = [r for r in all_rows if not r["excluded"]]
    
    if excluded:
        print(f"  Excluded from summary: {sorted(excluded)} (low coverage)")
    print(f"  Included in summary: {len(included_rows)} / {len(all_rows)} samples")
    
    # Write per-sample table
    per_sample_file = os.path.join(OUTPUT_DIR, f"bam_qc_per_sample_{output_prefix}.tsv")
    with open(per_sample_file, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["sample"] + columns + ["excluded"],
            delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for row in all_rows:
            out_row = {k: (v if v is not None else "NA") for k, v in row.items()}
            out_row["excluded"] = "yes" if row.get("excluded") else "no"
            writer.writerow(out_row)
    
    print(f"  Per-sample table: {per_sample_file}")
    print(f"  Total samples: {len(all_rows)}")
    
    # Compute and write summary
    summary = compute_summary_stats(included_rows, columns)
    summary_file = os.path.join(OUTPUT_DIR, f"bam_qc_summary_{output_prefix}.tsv")
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
    print(f"\n  === {output_prefix.upper()} BAM QC Summary ===")
    print(f"  {'Metric':<16} {'Mean':>10} {'SD':>8} {'Min':>10} {'Max':>10} {'Median':>10} {'N':>4}")
    print(f"  {'-'*16} {'-'*10} {'-'*8} {'-'*10} {'-'*10} {'-'*10} {'-'*4}")
    for col in columns:
        s = summary[col]
        vals = [
            f"{s['mean']}" if s['mean'] is not None else "NA",
            f"{s['sd']}" if s['sd'] is not None else "NA",
            f"{s['min']}" if s['min'] is not None else "NA",
            f"{s['max']}" if s['max'] is not None else "NA",
            f"{s['median']}" if s['median'] is not None else "NA",
        ]
        print(f"  {col:<16} {vals[0]:>10} {vals[1]:>8} {vals[2]:>10} {vals[3]:>10} {vals[4]:>10} {s['n']:>4}")


def main():
    print("=" * 60)
    print("BAM QC Table Generation (mosdepth + flagstat/samtools stats)")
    print("=" * 60)
    
    print(f"\nProcessing ILLUMINA GATK...")
    process_platform(ILLUMINA_GATK_DIR, "illumina_gatk", "illumina_gatk")
    
    print(f"\nProcessing AVITI GATK...")
    process_platform(AVITI_GATK_DIR, "aviti_gatk", "aviti_gatk")

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)


if __name__ == "__main__":
    main()
