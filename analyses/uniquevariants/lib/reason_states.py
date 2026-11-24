#!/usr/bin/env python3
"""reason_states.py - the genotype in the cohort raw GATK calls of each PPR-only allele that no raw call of the individual explains.
collect: the candidates of every individual and windows of 1 kb on both sides of them (targets.bed, Ensembl contig names).
states: reads bcftools query -H output on stdin; per individual and allele: carried, ref or missing (exact key after norm).
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uvlib  # noqa: E402

PAD = 1000
RANK = {"missing": 0, "ref": 1, "carried": 2}


def matched_keys(d):
    """the PPR keys matched (GM or AM) in any subset run of one eval dir"""
    lv = pd.read_csv(os.path.join(d, "subsets", "levels.tsv"), sep="\t", dtype=str)
    out = set()
    for name in lv.loc[lv["run"] == "yes", "level"]:
        r = pd.read_csv(os.path.join(d, f"classify_{name}.records.tsv.gz"), sep="\t", dtype=str, keep_default_na=False,
                        usecols=["side", "CHROM", "POS", "REF", "ALT", "class"])
        r = r[(r["side"] == "B") & r["class"].isin(["GM", "AM"])]
        out |= set(uvlib.make_key(r))
    return out


def inputs_of(d):
    return dict(ln.rstrip("\n").split("\t", 1) for ln in open(os.path.join(d, "inputs.tsv")) if "\t" in ln)


def ppr_unique(pair_records):
    r = pd.read_csv(pair_records, sep="\t", dtype=str, keep_default_na=False, usecols=["side", "CHROM", "POS", "REF", "ALT", "class"])
    r = r[(r["side"] == "B") & (r["class"] == "U")]
    return set(uvlib.make_key(r))


def collect(a):
    man = pd.read_csv(a.manifest, sep="\t", dtype=str, header=None, names=["patient", "dir"])
    rows = []
    for p, d in zip(man["patient"], man["dir"]):
        inp = inputs_of(d)
        u = ppr_unique(inp["pair_records"])
        raw = pd.read_csv(inp["raw_keys"], sep="\t", dtype=str, keep_default_na=False, usecols=["CHROM", "POS", "REF", "ALT"])
        cand = u - matched_keys(d) - set(uvlib.make_key(raw))
        rows += [(p, k) for k in sorted(cand)]
    df = pd.DataFrame(rows, columns=["patient", "key"])
    uvlib.write_tsv(df, a.out_prefix + ".candidates.tsv.gz")
    iv = {}
    for k in set(df["key"]):
        c, pos, _ = k.split(":", 2)
        c = c[3:] if c.startswith("chr") else c
        c = "MT" if c == "M" else c
        iv.setdefault(c, []).append((max(0, int(pos) - 1 - PAD), int(pos) + PAD))
    with open(a.out_prefix + ".targets.bed", "w") as f:
        for c in sorted(iv, key=lambda x: (len(x), x)):
            cur = None
            for s, e in sorted(iv[c]):
                if cur and s <= cur[1]:
                    cur[1] = max(cur[1], e)
                else:
                    if cur:
                        f.write(f"{c}\t{cur[0]}\t{cur[1]}\n")
                    cur = [s, e]
            if cur:
                f.write(f"{c}\t{cur[0]}\t{cur[1]}\n")
    print(f"reason_states collect: {len(df)} candidate (person, key) pairs, {df['key'].nunique()} distinct keys")


def states(a):
    cand = pd.read_csv(a.candidates, sep="\t", dtype=str, keep_default_na=False)
    want = {}
    for p, k in zip(cand["patient"], cand["key"]):
        want.setdefault(k, set()).add(p)
    lab = dict(ln.rstrip("\n").split("\t") for ln in open(a.label_map) if ln.strip())
    cols = None
    best = {}
    for line in sys.stdin:
        if line.startswith("#"):
            hdr = line.lstrip("# ").rstrip("\n").split("\t")
            # header fields look like "[2]<sample>:GT"
            cols = [lab.get(h.split("]", 1)[1].rsplit(":", 1)[0]) for h in hdr[1:]]
            continue
        f = line.rstrip("\n").split("\t")
        k = f[0]
        if k not in want:
            continue
        for i, gt in enumerate(f[1:]):
            p = cols[i]
            if p is None or p not in want[k]:
                continue
            al = gt.replace("|", "/").split("/")
            s = "carried" if any(x not in ("0", ".") for x in al) else ("missing" if "." in al else "ref")
            if RANK[s] > RANK.get(best.get((p, k)), -1):
                best[(p, k)] = s
    if cols is None:
        sys.exit("reason_states.py states: no header line on stdin (use bcftools query -H)")
    df = pd.DataFrame([(p, k, s) for (p, k), s in sorted(best.items())], columns=["patient", "key", "state"])
    uvlib.write_tsv(df, a.out)
    print(f"reason_states states: {len(df)} (person, key) pairs found in the raw merge, of {len(cand)} candidates")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["collect", "states"])
    ap.add_argument("--manifest")
    ap.add_argument("--out-prefix")
    ap.add_argument("--candidates")
    ap.add_argument("--label-map")
    ap.add_argument("--out")
    a = ap.parse_args()
    collect(a) if a.mode == "collect" else states(a)


if __name__ == "__main__":
    main()
