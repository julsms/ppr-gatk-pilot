#!/usr/bin/env python3
"""reason_subsets.py - the status of each raw GATK record of one individual in the final GATK file (PASS, VQSR, ExcessHet, OTHER,
ABSENT; ExcessHet wins over a VQSR tranche), and the cumulative subsets: vqsr = PASS + VQSR, excesshet = + ExcessHet,
other = + OTHER, absent = + ABSENT. A subset that adds no record is not written.
Outputs: <out-dir>/status.tsv.gz, levels.tsv, <level>.vcf.gz (+ .tbi), <level>.keys.tsv.gz, checks.tsv
"""
import argparse
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uvlib  # noqa: E402

LEVELS = [("vqsr", "a_vqsr", "VQSR"), ("excesshet", "a2_excesshet", "ExcessHet"), ("other", "a3_other_filter", "OTHER"),
          ("absent", "b_absent", "ABSENT")]
KEYCOLS = ["CHROM", "POS", "REF", "ALT", "GT", "TYPE", "SPANDEL", "AUTO"]


def status_of(flt):
    if flt == "PASS":
        return "PASS"
    vals = flt.split(";")
    if "ExcessHet" in vals:
        return "ExcessHet"
    if all(v.startswith("VQSRTranche") for v in vals):
        return "VQSR"
    return "OTHER"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-vcf", required=True)
    ap.add_argument("--raw-keys", required=True)
    ap.add_argument("--final-all-vcf", required=True)
    ap.add_argument("--final-pass-keys", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--bcftools", default=os.environ.get("BCFTOOLS", "bcftools"))
    ap.add_argument("--require-pass-equal", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    O = a.out_dir

    raw = pd.read_csv(a.raw_keys, sep="\t", dtype=str, keep_default_na=False, usecols=KEYCOLS)
    raw["key"] = uvlib.make_key(raw)
    fin = {}
    for line in uvlib.bcftools_lines(a.bcftools, ["query", "-f", "%CHROM:%POS:%REF:%ALT\t%FILTER\n", a.final_all_vcf]):
        k, f = line.rstrip("\n").split("\t")
        fin[k] = f
    fp = pd.read_csv(a.final_pass_keys, sep="\t", dtype=str, keep_default_na=False, usecols=["CHROM", "POS", "REF", "ALT"])
    fpk = set(uvlib.make_key(fp))
    raw["status"] = [status_of(fin[k]) if k in fin else "ABSENT" for k in raw["key"]]
    rawk = set(raw["key"])
    pass_raw = set(raw.loc[raw["status"] == "PASS", "key"])
    checks = [("raw_records", len(raw)), ("final_all_records", len(fin)), ("final_pass_records", len(fpk)),
              ("final_pass_not_in_raw", len(fpk - rawk)), ("final_all_not_in_raw", len(set(fin) - rawk)),
              ("raw_pass_not_in_final_pass", len(pass_raw - fpk)), ("final_pass_raw_status_not_pass", len(fpk - pass_raw)),
              ("require_pass_equal", int(a.require_pass_equal))]
    sc = raw["status"].value_counts()
    for s in ("PASS", "VQSR", "ExcessHet", "OTHER", "ABSENT"):
        checks.append((f"status_{s}", int(sc.get(s, 0))))

    # cumulative levels
    included = {"PASS"}
    levels, run_levels = [], {}
    for name, cls, st in LEVELS:
        included = included | {st}
        n_added = int(sc.get(st, 0))
        n_rec = int(raw["status"].isin(included).sum())
        run = n_added > 0
        levels.append({"level": name, "class": cls, "added_status": st, "n_added": n_added, "n_records": n_rec,
                       "run": "yes" if run else "no"})
        if run:
            run_levels[name] = set(included)
    st_of = dict(zip(raw["key"], raw["status"]))
    pipes = {}
    for name in run_levels:
        fo = open(os.path.join(O, f"{name}.vcf.gz.tmp"), "wb")
        pipes[name] = (fo, subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=fo, text=True))
    n_written = {n: 0 for n in run_levels}
    for line in uvlib.bcftools_lines(a.bcftools, ["view", a.raw_vcf]):
        if line.startswith("#"):
            for _, p in pipes.values():
                p.stdin.write(line)
            continue
        f = line.split("\t", 5)
        s = st_of[f"{f[0]}:{f[1]}:{f[3]}:{f[4]}"]
        for n, inc in run_levels.items():
            if s in inc:
                pipes[n][1].stdin.write(line)
                n_written[n] += 1
    for n, (fo, p) in pipes.items():
        p.stdin.close()
        if p.wait() != 0:
            sys.exit("bgzip failed")
        fo.close()
        os.replace(os.path.join(O, f"{n}.vcf.gz.tmp"), os.path.join(O, f"{n}.vcf.gz"))
        subprocess.run(["tabix", "-f", "-p", "vcf", os.path.join(O, f"{n}.vcf.gz")], check=True)
        sub = raw[raw["status"].isin(run_levels[n])][KEYCOLS]
        uvlib.write_tsv(sub, os.path.join(O, f"{n}.keys.tsv.gz"))
        lv = [x for x in levels if x["level"] == n][0]
        if n_written[n] != lv["n_records"]:
            sys.exit(f"reason_subsets.py: level {n}: {n_written[n]} records written, {lv['n_records']} expected")

    uvlib.write_tsv(raw[["key", "status"]], os.path.join(O, "status.tsv.gz"))
    pd.DataFrame(levels).to_csv(os.path.join(O, "levels.tsv"), sep="\t", index=False)
    with open(os.path.join(O, "checks.tsv"), "w") as f:
        f.write("check\tvalue\n")
        for c in checks:
            f.write(f"{c[0]}\t{c[1]}\n")
    if fpk - rawk:
        sys.exit(f"reason_subsets.py: {len(fpk - rawk)} final PASS records are not in the raw set; is the raw file the source?")
    if a.require_pass_equal and pass_raw - fpk:
        sys.exit(f"reason_subsets.py: {len(pass_raw - fpk)} raw records are PASS in the final file but not in the final PASS set "
                 "(the GATK set of the concordance); were the two prepared from the same file?")
    print("reason_subsets: " + ", ".join(f"{c[0]} {c[1]}" for c in checks if c[0].startswith("status_")))


if __name__ == "__main__":
    main()
