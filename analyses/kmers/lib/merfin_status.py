#!/usr/bin/env python3
"""merfin_status.py - Merfin status of every record of a call set (set: key, merfin 1 = kept by merfin -filter), and the counts of
one sample by region and type (sample: direct = merfin on the set itself; context = the unique records in the run of their full set).
Usage: python3 lib/merfin_status.py set --input X.vcf.gz --validated X.filter.vcf --out X.status.tsv.gz
       python3 lib/merfin_status.py sample --vcfs sheet/vcfs.tsv --platform illumina --sample <sample> --results <folder> --out <file>
"""
import argparse
import os
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "uniquevariants", "lib"))
import uvlib  # noqa: E402


def vcf_keys(bcftools, path):
    """the CHROM:POS:REF:ALT of every record (one ALT per record; a multi-ALT record would be an error)"""
    keys = []
    for line in uvlib.bcftools_lines(bcftools, ["query", "-f", "%CHROM:%POS:%REF:%ALT\n", path]):
        k = line.rstrip("\n")
        if "," in k.rsplit(":", 1)[1]:
            sys.exit(f"merfin_status.py: a record with more than one ALT allele in {path}; the inputs must be split")
        keys.append(k)
    return keys


def write_atomic(df, path):
    uvlib.write_tsv(df, path)


def mode_set(a):
    inp = vcf_keys(a.bcftools, a.input)
    if len(set(inp)) != len(inp):
        sys.exit(f"merfin_status.py: {len(inp) - len(set(inp))} duplicate keys in {a.input}")
    val = set(vcf_keys(a.bcftools, a.validated))
    extra = val - set(inp)
    if extra:
        sys.exit(f"merfin_status.py: {len(extra)} keys of the merfin output are not in the input {a.input}")
    df = pd.DataFrame({"key": inp})
    df["merfin"] = df["key"].isin(val).astype(int)
    write_atomic(df, a.out)
    print(f"{len(df)}\t{int(df['merfin'].sum())}")


def load_status(path):
    return pd.read_csv(path, sep="\t", dtype={"key": str, "merfin": int}, keep_default_na=False)


def done_fingerprint(path):
    """the fingerprint line of a DONE marker, or "" when the marker does not exist"""
    if not os.path.exists(path):
        return ""
    with open(path) as fh:
        for line in fh:
            k, _, v = line.rstrip("\n").partition("\t")
            if k == "fingerprint":
                return v
    return ""


EXIT_STATE = {"64": "missing_input", "65": "stale", "66": "no_disk"}


def sample_states(samples, results, state_path):
    """the state of every sample of the sheet: current when sample_state.tsv records exit 0 with the fingerprint of its SAMPLE_DONE"""
    st = {}
    if state_path and os.path.exists(state_path):
        t = pd.read_csv(state_path, sep="\t", dtype=str, keep_default_na=False)
        st = {(r.platform, r.sample): r for r in t.itertuples(index=False)}
    out = {}
    for r in samples.itertuples(index=False):
        fp = done_fingerprint(os.path.join(results, r.platform, r.sample, "SAMPLE_DONE"))
        s = st.get((r.platform, r.sample))
        if s is not None and s.exit == "0" and fp and fp == s.sample_done_fp:
            out[(r.platform, r.sample)] = "current"
        elif s is not None and s.exit != "0":
            out[(r.platform, r.sample)] = EXIT_STATE.get(s.exit, "failed")
        else:
            out[(r.platform, r.sample)] = "not_checked" if fp else "not_done"
    return out


def counts(df, strata, extra):
    c = uvlib.count_by(df, strata, extra=extra, value_cols=["merfin"]).rename(columns={"merfin": "validated"})
    c["rate"] = [round(v / n, 6) if n else float("nan") for v, n in zip(c["validated"], c["n"])]
    return c


def mode_sample(a):
    sheet = pd.read_csv(a.vcfs, sep="\t", dtype=str, keep_default_na=False)
    rows = sheet[(sheet["platform"] == a.platform) & (sheet["sample"] == a.sample)]
    if rows.empty:
        sys.exit(f"merfin_status.py: no row for {a.platform} {a.sample} in {a.vcfs}")
    strata = list(uvlib.strata_names())
    status = {}
    for r in rows.itertuples(index=False):
        if r.run == "yes" or r.kind == "all":
            p = os.path.join(a.results, f"{r.set_id}.status.tsv.gz")
            if not os.path.exists(p):
                sys.exit(f"merfin_status.py: missing {p}")
            status[r.set_id] = load_status(p)
    out = []
    for r in rows.itertuples(index=False):
        cols = ["CHROM", "POS", "REF", "ALT", "TYPE", "AUTO"]
        head = pd.read_csv(r.keys, sep="\t", nrows=0).columns
        keys = uvlib.load_keys(r.keys, usecols=cols + [s for s in strata if s in head])
        st = [s for s in strata if s in keys.columns]
        extra = {"platform": a.platform, "sample": a.sample, "set_id": r.set_id, "kind": r.kind, "pipeline": r.pipeline}
        if r.set_id in status:
            s = status[r.set_id]
            if len(s) != len(keys) or set(s["key"]) != set(keys["key"]):
                sys.exit(f"merfin_status.py: {r.set_id}: the key table {r.keys} and the merfin status hold different keys")
            m = keys.merge(s, on="key", how="left")
            out.append(counts(m, st, dict(extra, mode="direct")))
        if r.kind == "unique" and r.context not in ("", "."):
            if r.context not in status:
                sys.exit(f"merfin_status.py: {r.set_id}: no status of its full set {r.context}")
            c = status[r.context]
            m = keys.merge(c, on="key", how="left")
            lost = int(m["merfin"].isna().sum())
            if lost:
                sys.exit(f"merfin_status.py: {r.set_id}: {lost} unique keys are not in the full set {r.context}")
            m["merfin"] = m["merfin"].astype(int)
            out.append(counts(m, st, dict(extra, mode="context")))
    res = pd.concat(out, ignore_index=True)
    first = ["platform", "sample", "set_id", "kind", "pipeline", "mode", "stratum", "type", "n", "validated", "rate"]
    write_atomic(res[first], a.out)
    g = res[(res["stratum"] == "genome") & (res["type"] == "ALL")]
    for r in g.itertuples(index=False):
        print(f"{r.set_id}\t{r.mode}\t{r.validated}/{r.n}\t{100 * r.rate:.2f}%" if r.n else f"{r.set_id}\t{r.mode}\t0/0")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)
    s = sub.add_parser("set")
    s.add_argument("--input", required=True)
    s.add_argument("--validated", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--bcftools", default="bcftools")
    m = sub.add_parser("sample")
    m.add_argument("--vcfs", required=True)
    m.add_argument("--platform", required=True)
    m.add_argument("--sample", required=True)
    m.add_argument("--results", required=True)
    m.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.mode == "set":
        mode_set(a)
    else:
        mode_sample(a)


if __name__ == "__main__":
    main()
