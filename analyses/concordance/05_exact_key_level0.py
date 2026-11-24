#!/usr/bin/env python3
"""Exact-key concordance of one pair and individual on the same prepared files: records match on CHROM:POS:REF:ALT, genotypes on
the key and GT, with union denominators, per region and type (level0.tsv); and per record the key match against the vcfeval class
of 04_classify_count.py (bridge.tsv).
Usage: python3 05_exact_key_level0.py --pair <pair> --patient <individual>    (in the output folder, or with OUT set)
"""
import argparse
import csv
import datetime
import hashlib
import os
import resource
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "lib"))
import strata_names  # noqa: E402  (lib/strata_names.py: the stratum name rule and the column prefix)

TYPES = ["ALL", "SNV", "INDEL"]
EX_MISSING, EX_STALE = 64, 65


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def out_root():
    return os.environ.get("OUT", ".")


def pairs_table():
    return read_tsv(os.environ.get("PAIRS", os.path.join(HERE, "config", "pairs.tsv")))


def load_keys(path):
    """a key table of 01 with the 0/1 columns as booleans and every stratum column renamed stratum:<name>
    (lib/strata_names.py: no per-record column such as class can share a name with a stratum)"""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    strata = strata_names.check([c for c in df.columns if c not in ("CHROM", "POS", "REF", "ALT", "GT", "TYPE", "SPANDEL", "AUTO")],
                                f"05_exact_key_level0.py: {path}")
    df["key"] = df["CHROM"] + ":" + df["POS"] + ":" + df["REF"] + ":" + df["ALT"]
    df["gkey"] = df["key"] + ":" + df["GT"]
    for c in df.columns:
        if c not in ("CHROM", "POS", "REF", "ALT", "GT", "TYPE", "key", "gkey"):
            df[c] = df[c] == "1"
    return strata_names.namespaced(df, strata)


def strata_of(df):
    p = strata_names.COL_PREFIX
    return ["genome", "autosomes"] + [c[len(p):] for c in df.columns if c.startswith(p)]


def mask(df, st, t):
    m = pd.Series(True, index=df.index)
    if st == "autosomes":
        m &= df["AUTO"]
    elif st != "genome":
        m &= df[strata_names.COL_PREFIX + st]
    if t != "ALL":
        m &= df["TYPE"] == t
    return m


def level0(a, b, pair, patient):
    # a key has the same type and the same stratum membership in both files (same position, same alleles, same BED rule),
    # so the shared keys of a stratum are the keys of A in that stratum that B also has; union = N_A + N_B - shared
    in_b = a["key"].isin(pd.Index(b["key"]))
    g_in_b = a["gkey"].isin(pd.Index(b["gkey"]))
    rows = []
    for st in strata_of(a):
        for t in TYPES:
            ma, mb = mask(a, st, t), mask(b, st, t)
            na, nb = int(ma.sum()), int(mb.sum())
            shared, gts = int((ma & in_b).sum()), int((ma & g_in_b).sum())
            union = na + nb - shared
            rows.append({"patient": patient, "pair": pair, "stratum": st, "type": t, "L0_N_A": na, "L0_N_B": nb,
                         "L0_union": union, "L0_shared": shared, "L0_gt_shared": gts,
                         "L0_pos_conc": round(shared / union, 8) if union else float("nan"),
                         "L0_gt_conc": round(gts / union, 8) if union else float("nan"),
                         "L0_gt_rate": round(gts / shared, 8) if shared else float("nan"),
                         "L0_unique_A": na - shared, "L0_unique_B": nb - shared})
    return pd.DataFrame(rows)


def bridge(a, b, records, pair, patient):
    rec = pd.read_csv(records, sep="\t", dtype=str, keep_default_na=False)
    rec["key"] = rec["CHROM"] + ":" + rec["POS"] + ":" + rec["REF"] + ":" + rec["ALT"]
    vclass = {"GM": "shared", "AM": "shared", "U": "unique", "NA": "not_assessed"}
    rows = []
    for s, mine, theirs in (("A", a, b), ("B", b, a)):
        r = rec[rec["side"] == s][["key", "class"]]
        m = mine.merge(r, on="key", how="left")
        if m["class"].isna().any():
            sys.exit(f"05: side {s}: {int(m['class'].isna().sum())} prepared records without a class in {records}")
        m["key_match"] = m["key"].isin(set(theirs["key"])).map({True: "yes", False: "no"})
        m["vcfeval"] = m["class"].map(vclass)
        for st in strata_of(mine):
            for t in TYPES:
                sub = m[mask(m, st, t)]
                tab = sub.groupby(["key_match", "vcfeval"]).size()
                for km in ("yes", "no"):
                    for vc in ("shared", "unique", "not_assessed"):
                        rows.append({"patient": patient, "pair": pair, "side": s, "stratum": st, "type": t,
                                     "key_match": km, "vcfeval": vc, "records": int(tab.get((km, vc), 0))})
    return pd.DataFrame(rows)


def done_get(path, key):
    """a value of a DONE marker (key<TAB>value lines), or '.'"""
    if not path or not os.path.exists(path):
        return "."
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] == key:
                return p[1]
    return "."


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def one(a_keys, b_keys, records, prefix, pair, patient, fingerprint="", up=None):
    if os.path.exists(prefix + ".DONE"):
        if fingerprint and done_get(prefix + ".DONE", "fingerprint") == fingerprint:
            print(f"05 {pair} {patient}: already done (current)")
            sys.exit(0)
        print(f"STALE\tlevel0\t{pair}_{patient}\t{prefix}\ta prepared file, the classify table or the code changed")
        sys.exit(EX_STALE)
    for p in (a_keys, b_keys) + ((records,) if records else ()):
        if not os.path.exists(p):
            print(f"MISSING_INPUT\tlevel0\t{pair}\t{patient}\t{p}")
            sys.exit(EX_MISSING)
    os.makedirs(os.path.dirname(prefix), exist_ok=True)
    a, b = load_keys(a_keys), load_keys(b_keys)
    l0 = level0(a, b, pair, patient)
    l0.to_csv(prefix + ".level0.tsv.tmp", sep="\t", index=False, na_rep="NA")
    os.replace(prefix + ".level0.tsv.tmp", prefix + ".level0.tsv")
    if records:
        br = bridge(a, b, records, pair, patient)
        br.to_csv(prefix + ".bridge.tsv.tmp", sep="\t", index=False)
        os.replace(prefix + ".bridge.tsv.tmp", prefix + ".bridge.tsv")
    g = l0[(l0["stratum"] == "genome") & (l0["type"] == "ALL")].iloc[0]
    rss_mb = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss + resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss) / 1024
    done = [("date", datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")),
            ("fingerprint", fingerprint or "."), *(up or {}).items(),
            ("level0_md5", md5(prefix + ".level0.tsv")), ("peak_rss_mb", f"{rss_mb:.0f}"),
            ("genome_all", f"union {g['L0_union']} shared {g['L0_shared']} gt_shared {g['L0_gt_shared']}")]
    with open(prefix + ".DONE.tmp", "w") as f:
        f.write("".join(f"{k}\t{v}\n" for k, v in done))
    os.replace(prefix + ".DONE.tmp", prefix + ".DONE")
    print(f"05 {pair} {patient}: level 0 genome ALL N_A {g['L0_N_A']} N_B {g['L0_N_B']} shared {g['L0_shared']} "
          f"pos_conc {100 * g['L0_pos_conc']:.2f} % gt_conc {100 * g['L0_gt_conc']:.2f} %; peak RSS {rss_mb:.0f} MB")


def regression_compare(baseline_dir, expected_path):
    """--regression-compare: regression check of the development runs; its config files are not in this release"""
    root = out_root()
    smap = read_tsv(os.environ.get("SAMPLE_MAP", os.path.join(HERE, "config", "sample_map.tsv")))
    tabs = []
    for pr in pairs_table():
        d = os.path.join(root, "level0", pr["pair"])
        if not os.path.isdir(d):
            continue
        for f in sorted(os.listdir(d)):
            if f.endswith(".level0.tsv"):
                tabs.append(pd.read_csv(os.path.join(d, f), sep="\t"))
    if not tabs:
        sys.exit("05 --regression-compare: no level 0 tables under " + root)
    l0 = pd.concat(tabs)
    g = l0[l0["stratum"] == "genome"]
    out = []
    # 1. the expected values (--expected): means over patients and counts of one patient
    for e in read_tsv(expected_path):
        sub = g[(g["pair"] == e["pair"]) & (g["type"] == e["type"])]
        if e["kind"] == "mean_pct":
            val = 100 * sub[e["metric"]].mean() if len(sub) else float("nan")
            ok = abs(val - float(e["expected"])) <= float(e["tolerance"])
            out.append(("mean", e["pair"], ".", e["type"], e["metric"], f"{float(e['expected']):.4f}", f"{val:.4f}", len(sub),
                        "ok" if ok else "DIFFERENT"))
        elif e["kind"] == "count":
            sub = sub[sub["patient"] == e["patient"]]
            val = int(sub[e["metric"]].iloc[0]) if len(sub) else -1
            ok = val == int(e["expected"])
            out.append(("count", e["pair"], e["patient"], e["type"], e["metric"], e["expected"], val, 1, "ok" if ok else "DIFFERENT"))
    # 2. every sample of the per-platform tables of the development runs (--baseline-dir; PPR = side B, GATK = side A)
    if baseline_dir:
        tmap = {"SNPs": "SNV", "INDELs": "INDEL", "All": "ALL"}
        for plat in ("illumina", "aviti"):
            f = os.path.join(baseline_dir, f"PPR_GATK_{plat}", "final_summary", "all_samples_concordance.tsv")
            if not os.path.exists(f):
                out.append(("baseline", plat, ".", ".", "missing", f, ".", 0, "DIFFERENT"))
                continue
            base = pd.read_csv(f, sep="\t").drop_duplicates()
            col = f"{plat}_ppr"
            lab2pat = {r[col]: r["patient"] for r in smap}
            pair = f"GATK99_PPRold_{plat}"
            for _, r in base.iterrows():
                pat = lab2pat.get(str(r["Sample"]))
                sub = g[(g["pair"] == pair) & (g["patient"] == pat) & (g["type"] == tmap[r["Variant_Type"]])]
                for old, new in (("PPR_Variants", "L0_N_B"), ("GATK_Variants", "L0_N_A"), ("Position_Concordant", "L0_shared"),
                                 ("GT_Concordant", "L0_gt_shared")):
                    val = int(sub[new].iloc[0]) if len(sub) else -1
                    out.append(("sample", pair, pat, tmap[r["Variant_Type"]], new, int(r[old]), val, 1,
                                "ok" if val == int(r[old]) else "DIFFERENT"))
    res = pd.DataFrame(out, columns=["check", "pair", "patient", "type", "metric", "expected", "observed", "n", "verdict"])
    res.to_csv(os.path.join(root, "regression_comparison.tsv"), sep="\t", index=False)
    bad = res[res["verdict"] != "ok"]
    print(f"05 --regression-compare: {len(res)} checks, {len(bad)} different; table {root}/regression_comparison.tsv")
    print(res[res["check"] != "sample"].to_string(index=False))
    sys.exit(1 if len(bad) else 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair")
    ap.add_argument("--patient")
    ap.add_argument("--a-keys")
    ap.add_argument("--b-keys")
    ap.add_argument("--records")
    ap.add_argument("--out-prefix")
    ap.add_argument("--no-bridge", action="store_true")
    ap.add_argument("--regression-compare", action="store_true")
    ap.add_argument("--baseline-dir")
    ap.add_argument("--expected", default=os.path.join(HERE, "config", "regression_expected.tsv"))
    ap.add_argument("--fingerprint", default="")
    a = ap.parse_args()

    if a.regression_compare:
        regression_compare(a.baseline_dir, a.expected)
        return
    up = {}
    if a.pair and a.patient and not a.a_keys:
        r = [x for x in pairs_table() if x["pair"] == a.pair]
        if not r:
            sys.exit(f"05: unknown pair {a.pair}")
        r = r[0]
        root = out_root()
        a_dir = os.path.join(root, "prep", f"{r['a_set']}_{r['a_platform']}")
        b_dir = os.path.join(root, "prep", f"{r['b_set']}_{r['b_platform']}")
        a_keys = os.path.join(a_dir, f"{a.patient}.keys.tsv.gz")
        b_keys = os.path.join(b_dir, f"{a.patient}.keys.tsv.gz")
        records = None if a.no_bridge else os.path.join(root, "classify", a.pair, f"{a.patient}.records.tsv.gz")
        prefix = os.path.join(root, "level0", a.pair, a.patient)
        # the fingerprints of the upstream DONE markers (06 checks this chain)
        up = {"a_prep_fp": done_get(os.path.join(a_dir, f"{a.patient}.DONE"), "fingerprint"),
              "b_prep_fp": done_get(os.path.join(b_dir, f"{a.patient}.DONE"), "fingerprint"),
              "classify_fp": "nobridge" if a.no_bridge else
              done_get(os.path.join(root, "classify", a.pair, f"{a.patient}.DONE"), "fingerprint")}
    else:
        a_keys, b_keys, records, prefix = a.a_keys, a.b_keys, (None if a.no_bridge else a.records), a.out_prefix
    one(a_keys, b_keys, records, prefix, a.pair or ".", a.patient or ".", a.fingerprint, up)


if __name__ == "__main__":
    main()
