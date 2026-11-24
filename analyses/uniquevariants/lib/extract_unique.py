#!/usr/bin/env python3
"""extract_unique.py - the unique set of one individual: the records of one side of a pair with class U in the vcfeval classes
of the concordance step (for two pairs, U in both). The records are taken from the prepared VCF by key; the count per region and
type must equal U of that side in the concordance counts.
Outputs: <prefix>.vcf.gz (+ .tbi), <prefix>.keys.tsv.gz, <prefix>.counts.tsv, <prefix>.checks.tsv
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uvlib  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--side", required=True, choices=["A", "B"])
    ap.add_argument("--keys", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--uset", default=".")
    ap.add_argument("--patient", default=".")
    ap.add_argument("--bcftools", default=os.environ.get("BCFTOOLS", "bcftools"))
    a = ap.parse_args()
    recs, cnts = a.records.split(","), a.counts.split(",")
    if len(recs) != len(cnts):
        sys.exit("extract_unique.py: --records and --counts differ in length")

    fail, checks = [], []
    u_sets, lm_sets = [], []
    for r in recs:
        df = pd.read_csv(r, sep="\t", dtype=str, keep_default_na=False,
                         usecols=["side", "CHROM", "POS", "REF", "ALT", "class", "BK"])
        df = df[df["side"] == a.side]
        df["key"] = uvlib.make_key(df)
        u = df[df["class"] == "U"]
        u_sets.append(set(u["key"]))
        lm_sets.append(set(u.loc[u["BK"] == "lm", "key"]))
        checks.append((os.path.basename(os.path.dirname(r)), "U_in_records", len(u_sets[-1])))
    uk = set.intersection(*u_sets)
    lm = set.union(*lm_sets) & uk

    keys = uvlib.load_keys(a.keys)
    sel = keys[keys["key"].isin(uk)].copy()
    if len(sel) != len(uk):
        fail.append(f"{len(uk) - len(sel)} U keys are not in the prepared key table {a.keys}")
    sel["BK_lm"] = sel["key"].isin(lm)          # the lm flag of 04 (LMFLAG_A / LMFLAG_B there); LM is the low-mappability stratum
    strata = uvlib.strata_columns(keys)
    cnt = uvlib.count_by(sel, strata, extra={"uset": a.uset, "patient": a.patient}, value_cols=["BK_lm"])

    if len(recs) == 1:
        conc = pd.read_csv(cnts[0], sep="\t", dtype=str, keep_default_na=False)
        m = cnt.merge(conc[["stratum", "type", f"U_{a.side}", f"LMFLAG_{a.side}"]], on=["stratum", "type"], how="left")
        checks.append((".", "strata_not_compared", ",".join(uvlib.CONC_SKIP_STRATA) or "none"))
        m = m[~m["stratum"].isin(uvlib.CONC_SKIP_STRATA)]
        bad = m[(m[f"U_{a.side}"].astype(int) != m["n"]) | (m[f"LMFLAG_{a.side}"].astype(int) != m["BK_lm"])]
        checks.append((".", "rows_equal_to_concordance_U", len(m) - len(bad)))
        checks.append((".", "rows_differing_from_concordance_U", len(bad)))
        if len(bad):
            fail.append(f"{len(bad)} stratum x type rows differ from U/LMFLAG of {cnts[0]}")
    g = cnt[(cnt["stratum"] == "genome") & (cnt["type"] == "ALL")].iloc[0]

    os.makedirs(os.path.dirname(a.out_prefix), exist_ok=True)
    n_vcf = uvlib.filter_vcf_by_keys(a.bcftools, a.vcf, set(sel["key"]), a.out_prefix + ".vcf.gz")
    if n_vcf != len(sel):
        fail.append(f"the VCF has {n_vcf} records, the key table {len(sel)}")
    checks += [(".", "unique_records", len(sel)), (".", "lm_records", int(sel["BK_lm"].sum())), (".", "vcf_records", n_vcf)]

    out = sel.drop(columns=["key"])
    for c in out.columns:
        if out[c].dtype == bool:
            out[c] = out[c].astype(int)
    uvlib.write_tsv(out, a.out_prefix + ".keys.tsv.gz")
    uvlib.write_tsv(cnt, a.out_prefix + ".counts.tsv")
    with open(a.out_prefix + ".checks.tsv", "w") as f:
        f.write("source\tcheck\tvalue\n")
        for c in checks:
            f.write("\t".join(str(x) for x in c) + "\n")
        for x in fail:
            f.write(f".\tFAIL\t{x}\n")
    if fail:
        sys.exit("extract_unique.py: " + "; ".join(fail))
    print(f"extract_unique {a.uset} {a.patient}: {int(g['n'])} unique records ({int(g['BK_lm'])} lm)")


if __name__ == "__main__":
    main()
