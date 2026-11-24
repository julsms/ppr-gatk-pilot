#!/usr/bin/env python3
"""reason_classes.py - the GATK status of every PPR-only allele (vs GATK) of one individual: the first class that applies of
a_vqsr, a2_excesshet, a3_other_filter, b_absent (vcfeval matches it in that subset of the raw calls), a_unmatched (a raw call
with the same key), c_ref, c_missing (the allele is in the cohort raw calls, this individual is 0/0 or missing), d_none.
The classes must add up to U of the PPR side in every region and type.
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reason_states  # noqa: E402
import uvlib  # noqa: E402

CLASSES = ["a_vqsr", "a2_excesshet", "a3_other_filter", "b_absent", "a_unmatched", "c_ref", "c_missing", "d_none"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eval-dir", required=True)
    ap.add_argument("--states", required=True)
    ap.add_argument("--patient", required=True)
    ap.add_argument("--mode", default=".")
    ap.add_argument("--platform", default=".")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    d = a.eval_dir
    inp = reason_states.inputs_of(d)

    u = reason_states.ppr_unique(inp["pair_records"])
    keys = uvlib.load_keys(inp["ppr_keys"])
    sel = keys[keys["key"].isin(u)].copy()
    if len(sel) != len(u):
        sys.exit(f"reason_classes.py: {len(u) - len(sel)} PPR-unique keys are not in the PPR key table")
    cls = pd.Series("", index=sel["key"].values)

    lv = pd.read_csv(os.path.join(d, "subsets", "levels.tsv"), sep="\t", dtype=str)
    for name, c in zip(lv.loc[lv["run"] == "yes", "level"], lv.loc[lv["run"] == "yes", "class"]):
        r = pd.read_csv(os.path.join(d, f"classify_{name}.records.tsv.gz"), sep="\t", dtype=str, keep_default_na=False,
                        usecols=["side", "CHROM", "POS", "REF", "ALT", "class"])
        r = r[(r["side"] == "B") & r["class"].isin(["GM", "AM"])]
        m = set(uvlib.make_key(r))
        hit = (cls == "") & cls.index.isin(m)
        cls[hit] = c
    raw = pd.read_csv(inp["raw_keys"], sep="\t", dtype=str, keep_default_na=False, usecols=["CHROM", "POS", "REF", "ALT"])
    hit = (cls == "") & cls.index.isin(set(uvlib.make_key(raw)))
    cls[hit] = "a_unmatched"
    st = pd.read_csv(a.states, sep="\t", dtype=str, keep_default_na=False)
    st = st[st["patient"] == a.patient].set_index("key")["state"]
    s = st.reindex(cls.index)
    cls[(cls == "") & (s == "carried")] = "a_unmatched"
    cls[(cls == "") & (s == "ref")] = "c_ref"
    cls[(cls == "") & (s == "missing")] = "c_missing"
    cls[cls == ""] = "d_none"
    sel["class"] = cls.reindex(sel["key"]).values

    strata = uvlib.strata_columns(keys)
    cnt = uvlib.class_counts(sel, strata, "class", CLASSES,
                             extra={"mode": a.mode, "platform": a.platform, "patient": a.patient})
    conc = pd.read_csv(inp["pair_counts"], sep="\t", dtype=str, keep_default_na=False)
    m = cnt.merge(conc[["stratum", "type", "U_B"]], on=["stratum", "type"], how="left")
    m = m[~m["stratum"].isin(uvlib.CONC_SKIP_STRATA)]          # see uvlib.CONC_SKIP_STRATA
    bad = m[m["U_B"].astype(int) != m["n"]]
    fail = [f"{len(bad)} stratum x type rows: the classes do not add up to U_B of {inp['pair_counts']}"] if len(bad) else []
    out = sel[["CHROM", "POS", "REF", "ALT", "TYPE", "class"] + strata + ["AUTO"]].copy()
    for c in out.columns:
        if out[c].dtype == bool:
            out[c] = out[c].astype(int)
    uvlib.write_tsv(out, a.out_prefix + ".classes.tsv.gz")
    uvlib.write_tsv(cnt, a.out_prefix + ".counts.tsv")
    g = cnt[(cnt["stratum"] == "genome") & (cnt["type"] == "ALL")].iloc[0]
    with open(a.out_prefix + ".checks.tsv", "w") as f:
        f.write("check\tvalue\n")
        f.write(f"ppr_unique\t{len(sel)}\nrows_equal_to_U_B\t{len(m) - len(bad)}\n")
        for c in CLASSES:
            f.write(f"genome_{c}\t{int(g[c])}\n")
        for x in fail:
            f.write(f"FAIL\t{x}\n")
    if fail:
        sys.exit("reason_classes.py: " + "; ".join(fail))
    print(f"reason_classes {a.mode} {a.platform} {a.patient}: " + ", ".join(f"{c} {int(g[c])}" for c in CLASSES))


if __name__ == "__main__":
    main()
