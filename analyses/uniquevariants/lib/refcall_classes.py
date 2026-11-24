#!/usr/bin/env python3
"""refcall_classes.py - the RefCall class of each unique allele against a DeepVariant call set: i_refcall when its key is a
RefCall record (a candidate that DeepVariant rejected) of the other call set of the individual, ii_none otherwise.
Inputs: the key table of the unique set and the RefCall keys (CHROM POS REF ALT, no header).
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uvlib  # noqa: E402

CLASSES = ["i_refcall", "ii_none"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--unique-keys", required=True)
    ap.add_argument("--refcall-keys", required=True)
    ap.add_argument("--uset", default=".")
    ap.add_argument("--patient", default=".")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    u = uvlib.load_keys(a.unique_keys)
    rc = set()
    with open(a.refcall_keys) as f:
        for ln in f:
            c, p, r, x = ln.rstrip("\n").split("\t")
            rc.add(f"{c}:{p}:{r}:{x}")
    u["class"] = ["i_refcall" if k in rc else "ii_none" for k in u["key"]]
    strata = uvlib.strata_columns(u)
    cnt = uvlib.class_counts(u, strata, "class", CLASSES, extra={"uset": a.uset, "patient": a.patient})
    uc = uvlib.count_by(u, strata)
    ok = (cnt["n"].values == uc["n"].values).all()
    out = u[["CHROM", "POS", "REF", "ALT", "TYPE", "class"]]
    uvlib.write_tsv(out, a.out_prefix + ".classes.tsv.gz")
    uvlib.write_tsv(cnt, a.out_prefix + ".counts.tsv")
    g = cnt[(cnt["stratum"] == "genome") & (cnt["type"] == "ALL")].iloc[0]
    with open(a.out_prefix + ".checks.tsv", "w") as f:
        f.write("check\tvalue\n")
        f.write(f"unique\t{len(u)}\nrefcall_keys\t{len(rc)}\ni_refcall\t{int(g['i_refcall'])}\nii_none\t{int(g['ii_none'])}\n")
        f.write(f"closure\t{'ok' if ok else 'FAIL'}\n")
    if not ok:
        sys.exit("refcall_classes.py: the classes do not add up to the unique count")
    print(f"refcall_classes {a.uset} {a.patient}: i_refcall {int(g['i_refcall'])}, ii_none {int(g['ii_none'])}")


if __name__ == "__main__":
    main()
