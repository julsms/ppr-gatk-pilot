#!/usr/bin/env python3
"""BAM QC of the PPR alignments (projected onto GRCh38) per sample and per platform, as for GATK: mean coverage, % of the genome
at 10x, 20x and 30x, % mapped and % properly paired reads.
Input: bam_qc_ppr/{illumina_ppr,aviti_ppr}/ (bam_qc.sh: <sample>.mosdepth.summary.txt, .mosdepth.global.dist.txt, _flagstat.txt).
Output: tables_PPRbam_qc/bam_qc_per_sample_<platform>_ppr.tsv and bam_qc_summary_<platform>_ppr.tsv.
"""

import os
import re
import glob
import csv
import statistics

# ============================================================
# CONFIGURATION
# ============================================================
PPR_BASE = "bam_qc_ppr"
ILLUMINA_PPR_DIR = os.path.join(PPR_BASE, "illumina_ppr")
AVITI_PPR_DIR = os.path.join(PPR_BASE, "aviti_ppr")
OUTPUT_DIR = "tables_PPRbam_qc"

EXCLUDED_SAMPLES = {
    "illumina_ppr": ["23"],
    "aviti_ppr": ["23"],
}

# ============================================================
# PARSING FUNCTIONS
# ============================================================

def parse_mosdepth_summary(filepath):
    """Get mean coverage from mosdepth summary (total row)."""
    if not os.path.isfile(filepath):
        return None
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4 and parts[0] in ("total", "total_region"):
                try:
                    return float(parts[3])
                except (ValueError, IndexError):
                    pass
    return None


def parse_mosdepth_global_dist(filepath):
    """Get % genome >=10x, >=20x, >=30x from global.dist.txt."""
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


def parse_flagstat(filepath):
    """Get % mapped from flagstat."""
    if not os.path.isfile(filepath):
        return None
    with open(filepath, "r") as f:
        content = f.read()
    m = re.search(r"(\d+)\s+\+\s+\d+\s+mapped\s+\((\d+\.?\d*)%", content)
    if m:
        return float(m.group(2))
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

# ============================================================
# SAMPLE DISCOVERY AND QC EXTRACTION
# ============================================================

def discover_samples(qc_dir):
    """Find all samples in a flat mosdepth output directory.
    
    Files are named: SAMPLE.mosdepth.summary.txt, SAMPLE.mosdepth.global.dist.txt, SAMPLE_flagstat.txt
    """
    samples = set()
    for f in os.listdir(qc_dir):
        if f.endswith(".mosdepth.summary.txt"):
            sample = f.replace(".mosdepth.summary.txt", "")
            samples.add(sample)
    return sorted(samples)


def extract_bam_qc(qc_dir, sample_name):
    """Extract BAM QC for one sample from flat dir."""
    row = {"sample": sample_name}
    
    summary_file = os.path.join(qc_dir, f"{sample_name}.mosdepth.summary.txt")
    dist_file = os.path.join(qc_dir, f"{sample_name}.mosdepth.global.dist.txt")
    flagstat_file = os.path.join(qc_dir, f"{sample_name}_flagstat.txt")
    
    row["mean_coverage"] = parse_mosdepth_summary(summary_file)
    
    thresholds = parse_mosdepth_global_dist(dist_file)
    row["pct_10x"] = thresholds["pct_10x"]
    row["pct_20x"] = thresholds["pct_20x"]
    row["pct_30x"] = thresholds["pct_30x"]
    
    row["pct_mapped"] = parse_flagstat(flagstat_file)
    row["pct_properly_paired"] = parse_flagstat_properly_paired(flagstat_file)
    
    return row


# ============================================================
# MAIN
# ============================================================

def compute_summary_stats(rows, columns):
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


def process_platform(qc_dir, output_prefix, platform_key):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    if not os.path.isdir(qc_dir):
        print(f"  ERROR: Directory not found: {qc_dir}")
        return
    
    samples = discover_samples(qc_dir)
    if not samples:
        print(f"  WARNING: No samples found in {qc_dir}")
        return
    
    columns = ["mean_coverage", "pct_10x", "pct_20x", "pct_30x", "pct_mapped", "pct_properly_paired"]
    excluded = set(EXCLUDED_SAMPLES.get(platform_key, []))
    
    all_rows = []
    for sample in samples:
        row = extract_bam_qc(qc_dir, sample)
        row["excluded"] = sample in excluded
        all_rows.append(row)
    
    included_rows = [r for r in all_rows if not r["excluded"]]
    
    print(f"  Samples found: {len(all_rows)}")
    if excluded:
        print(f"  Excluded from summary: {sorted(excluded)}")
    print(f"  Included in summary: {len(included_rows)}")
    
    # Per-sample TSV
    per_sample_file = os.path.join(OUTPUT_DIR, f"bam_qc_per_sample_{output_prefix}.tsv")
    with open(per_sample_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample"] + columns + ["excluded"],
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in all_rows:
            out_row = {k: (v if v is not None else "NA") for k, v in row.items()}
            out_row["excluded"] = "yes" if row.get("excluded") else "no"
            writer.writerow(out_row)
    print(f"  Per-sample: {per_sample_file}")
    
    # Summary TSV
    summary = compute_summary_stats(included_rows, columns)
    summary_file = os.path.join(OUTPUT_DIR, f"bam_qc_summary_{output_prefix}.tsv")
    with open(summary_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["metric", "mean", "sd", "min", "max", "median", "n"])
        for col in columns:
            s = summary[col]
            writer.writerow([col, s["mean"], s["sd"], s["min"], s["max"], s["median"], s["n"]])
    print(f"  Summary: {summary_file}")
    
    # Print
    print(f"\n  === {output_prefix.upper()} BAM QC Summary ===")
    print(f"  {'Metric':<16} {'Mean':>10} {'SD':>8} {'Min':>10} {'Max':>10} {'Median':>10}")
    print(f"  {'-'*70}")
    for col in columns:
        s = summary[col]
        if s["mean"] is not None:
            print(f"  {col:<16} {s['mean']:>10} {s['sd']:>8} {s['min']:>10} {s['max']:>10} {s['median']:>10}")
        else:
            print(f"  {col:<16} {'NA':>10}")


def main():
    print("=" * 60)
    print("BAM QC Tables - PPR Illumina + PPR Aviti")
    print("=" * 60)
    
    print(f"\nProcessing ILLUMINA PPR...")
    process_platform(ILLUMINA_PPR_DIR, "illumina_ppr", "illumina_ppr")
    
    print(f"\nProcessing AVITI PPR...")
    process_platform(AVITI_PPR_DIR, "aviti_ppr", "aviti_ppr")
    
    print("\nDONE")


if __name__ == "__main__":
    main()
