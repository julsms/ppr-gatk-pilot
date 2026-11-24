#!/usr/bin/env python3
"""Cohort tables of the concordance from the classify/ and level0/ tables of 04 and 05: all_counts.tsv (one row per individual,
pair, region and type), genome_counts.tsv (genome rows), pair_means.tsv, checks.tsv, and paired Wilcoxon tests (wilcoxon_H1-H3.tsv).
Usage: python3 06_summary_stats.py    (in the output folder, or with OUT set)
"""
import csv
import glob
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get("OUT", ".")
S = os.path.join(OUT, "summary")
METRICS = ["J_gt", "J_al"]
TEST_TYPES = ["ALL", "SNV", "INDEL"]


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(df, name):
    os.makedirs(S, exist_ok=True)
    p = os.path.join(S, name)
    df.to_csv(p + ".tmp", sep="\t", index=False, na_rep="NA")
    os.replace(p + ".tmp", p)


def collect(sub, suffix):
    tabs = []
    for f in sorted(glob.glob(os.path.join(OUT, sub, "*", f"*{suffix}"))):
        if os.path.exists(f[: -len(suffix)] + ".DONE"):
            tabs.append(pd.read_csv(f, sep="\t", dtype={"patient": str}))
    return pd.concat(tabs, ignore_index=True) if tabs else pd.DataFrame()


def done_get(path, key):
    """a value of a DONE marker (key<TAB>value lines), or '.'"""
    if not os.path.exists(path):
        return "."
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] == key:
                return p[1]
    return "."


def chain_check(pairs):
    """every classify and level 0 DONE must name the fingerprints of the upstream DONE markers that exist now"""
    bad = []
    for pr in pairs:
        a_dir = os.path.join(OUT, "prep", f"{pr['a_set']}_{pr['a_platform']}")
        b_dir = os.path.join(OUT, "prep", f"{pr['b_set']}_{pr['b_platform']}")
        up_dir = {"classify": ("vcfeval_fp", os.path.join(OUT, "vcfeval", pr["pair"])),
                  "level0": ("classify_fp", os.path.join(OUT, "classify", pr["pair"]))}
        for stage, (link, link_dir) in up_dir.items():
            for d in sorted(glob.glob(os.path.join(OUT, stage, pr["pair"], "*.DONE"))):
                pat = os.path.basename(d)[: -len(".DONE")]
                need = [("a_prep_fp", os.path.join(a_dir, pat + ".DONE")), ("b_prep_fp", os.path.join(b_dir, pat + ".DONE")),
                        (link, os.path.join(link_dir, pat + ".DONE"))]
                for key, up in need:
                    rec = done_get(d, key)
                    if key == "classify_fp" and rec == "nobridge":
                        continue
                    if rec != done_get(up, "fingerprint"):
                        bad.append(f"{stage}/{pr['pair']}/{pat}:{key}")
    return len(bad), bad


def ledger_value(s, p, patient, key):
    f = os.path.join(OUT, "prep", f"{s}_{p}", f"{patient}.ledger.tsv")
    if not os.path.exists(f):
        return np.nan
    for r in read_tsv(f):
        if r["step"] == key:
            return int(r["count"])
    return np.nan


def wilcoxon(x, y):
    d = np.asarray(x) - np.asarray(y)
    d = d[~np.isnan(d)]
    n_all = len(d)
    nz = d[d != 0]
    if len(nz) < 3:
        return n_all, np.nan, np.nan, np.nan
    ranks = stats.rankdata(np.abs(nz))
    v = float(ranks[nz > 0].sum())
    res = stats.wilcoxon(nz, zero_method="wilcox", correction=True, method="approx")
    p = float(res.pvalue)
    z = stats.norm.isf(p / 2) if p > 0 else np.inf
    return n_all, v, p, min(abs(z) / math.sqrt(n_all), 1.0)  # r = |z| / sqrt(n), n = all patients


def holm(p):
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    idx = np.where(ok)[0]
    if len(idx) == 0:
        return out
    order = idx[np.argsort(p[idx])]
    m = len(order)
    run = 0.0
    for i, k in enumerate(order):
        run = max(run, min(1.0, (m - i) * p[k]))
        out[k] = run
    return out


def tests(wide, family, comparisons):
    """comparisons: list of (label, x-spec, y-spec); a spec is a list of pairs whose mean is taken per patient"""
    rows = []
    for t in TEST_TYPES:
        for metric in METRICS:
            block = []
            for st in sorted(wide["stratum"].unique()):
                w = wide[(wide["stratum"] == st) & (wide["type"] == t)]
                for label, xs, ys in comparisons:
                    cols = [f"{metric}|{p}" for p in xs + ys]
                    if not all(c in w.columns for c in cols):
                        continue
                    ww = w.dropna(subset=cols)
                    x = ww[[f"{metric}|{p}" for p in xs]].mean(axis=1)
                    y = ww[[f"{metric}|{p}" for p in ys]].mean(axis=1)
                    n, v, p, r = wilcoxon(x, y)
                    block.append({"family": family, "comparison": label, "stratum": st, "type": t, "metric": metric, "n": n,
                                  "mean_x": x.mean(), "mean_y": y.mean(), "mean_diff": (x - y).mean(), "V": v, "p_raw": p,
                                  "effect_r": r})
            ph = holm([b["p_raw"] for b in block])
            for b, q in zip(block, ph):
                b["p_holm"] = q
                b["significant"] = bool(q < 0.05) if not np.isnan(q) else False
            rows += block
    return pd.DataFrame(rows)


def main():
    pairs = read_tsv(os.environ.get("PAIRS", os.path.join(HERE, "config", "pairs.tsv")))
    smap = read_tsv(os.environ.get("SAMPLE_MAP", os.path.join(HERE, "config", "sample_map.tsv")))
    patients = [r["patient"] for r in smap if r["exclude"] != "yes"]
    cnt = collect("classify", ".counts.tsv")
    l0 = collect("level0", ".level0.tsv")
    if cnt.empty:
        sys.exit("06: no finished classify tables under " + OUT)
    keys = ["patient", "pair", "stratum", "type"]
    allc = cnt.merge(l0, on=keys, how="left") if not l0.empty else cnt
    # pre-split PASS records (ledger of 01; NA without a ledger), per side
    pre = []
    for pr in pairs:
        for pat in allc.loc[allc["pair"] == pr["pair"], "patient"].unique():
            pre.append({"patient": pat, "pair": pr["pair"],
                        "presplit_PASS_A": ledger_value(pr["a_set"], pr["a_platform"], pat, "pass_records_presplit"),
                        "presplit_PASS_B": ledger_value(pr["b_set"], pr["b_platform"], pat, "pass_records_presplit")})
    allc = allc.merge(pd.DataFrame(pre), on=["patient", "pair"], how="left")
    allc["platform"] = [r["a_platform"] if r["kind"] == "within" else "cross" for _, r in allc.iterrows()]
    write(allc, "all_counts.tsv")

    r25_cols = ["patient", "platform", "pair", "type", "N_A", "N_B", "S_A", "S_B", "GM_A", "GM_B", "AM_A", "AM_B", "U_A", "U_B",
                "LMFLAG_A", "LMFLAG_B", "NA_A", "NA_B", "presplit_PASS_A", "presplit_PASS_B", "J_al", "J_gt", "share_A", "share_B",
                "gt_agree", "rep_gap", "NA_share", "J_gt_nospandel", "L0_union", "L0_shared", "L0_gt_shared", "L0_pos_conc",
                "L0_gt_conc", "L0_unique_A", "L0_unique_B"]
    g = allc[allc["stratum"] == "genome"]
    write(g[[c for c in r25_cols if c in g.columns]].sort_values(["pair", "type", "patient"]), "genome_counts.tsv")

    mcols = [c for c in ["J_al", "J_gt", "share_A", "share_B", "gt_agree", "J_gt_nospandel", "L0_pos_conc", "L0_gt_conc"]
             if c in allc.columns]
    agg = allc.groupby(["pair", "stratum", "type"])[mcols].agg(["count", "mean", "std", "min", "max"])
    agg.columns = [f"{a}_{b}" for a, b in agg.columns]
    agg = agg.reset_index()
    agg["NA_share_max"] = allc.groupby(["pair", "stratum", "type"])["NA_share"].max().values
    write(agg, "pair_means.tsv")

    # ---- checks ----
    chk, fail = [], []
    for pr in pairs:
        have = sorted(g.loc[g["pair"] == pr["pair"], "patient"].unique())
        avail = [p for p in patients if os.path.exists(os.path.join(OUT, "prep", f"{pr['a_set']}_{pr['a_platform']}", f"{p}.DONE"))
                 and os.path.exists(os.path.join(OUT, "prep", f"{pr['b_set']}_{pr['b_platform']}", f"{p}.DONE"))]
        gate_dir = os.path.join(OUT, "gate", pr["pair"])
        gate_ok = os.environ.get("SKIP_GATE") == "1" or os.path.exists(os.path.join(gate_dir, "GATE_PASSED"))
        note = "."
        if not gate_ok:
            avail, note = [], "gate not passed (waiting or failed): no rows expected"
        missing = sorted(set(avail) - set(have))
        chk.append({"check": "rows_per_pair", "pair": pr["pair"], "value": len(have),
                    "expected": len(avail), "of_cohort": len(patients),
                    "detail": ("missing: " + ",".join(missing)) if missing else note})
        if missing:
            fail.append(f"{pr['pair']}: {len(missing)} patients prepared but not counted")
    for s in ("A", "B"):
        bad = allc[allc[f"N_{s}"] != allc[f"GM_{s}"] + allc[f"AM_{s}"] + allc[f"U_{s}"] + allc[f"NA_{s}"]]
        chk.append({"check": f"closure_{s}", "pair": "all", "value": len(bad), "expected": 0, "of_cohort": ".", "detail": "."})
        if len(bad):
            fail.append(f"closure fails in {len(bad)} rows (side {s})")
    sbar = (allc["S_A"] + allc["S_B"]) / 2
    back = sbar / allc["J_al"] - sbar
    badb = ((back - (allc["U_A"] + allc["U_B"])).abs() > 1e-4 * (sbar + allc["U_A"] + allc["U_B"] + 1)) & allc["J_al"].gt(0)
    chk.append({"check": "back_calculation_U", "pair": "all", "value": int(badb.sum()), "expected": 0, "of_cohort": ".",
                "detail": "U_A + U_B = Sbar / J_al - Sbar"})
    if badb.any():
        fail.append("back-calculation of U fails")
    ga = g[g["type"] == "ALL"]
    alarms = ga[ga["NA_share"] > 0.001]
    chk.append({"check": "NA_alarm_genome", "pair": "all", "value": len(alarms), "expected": 0, "of_cohort": ".",
                "detail": ";".join(f"{r.pair}/{r.patient}:{r.NA_share:.5f}" for r in alarms.itertuples()) or "."})
    nb, detail = chain_check(pairs)
    chk.append({"check": "done_chain", "pair": "all", "value": nb, "expected": 0, "of_cohort": ".",
                "detail": ";".join(detail[:20]) or "every classify and level 0 table was made from the current upstream files"})
    if nb:
        fail.append(f"{nb} classify or level 0 tables were made from upstream files that changed since (stale; see checks.tsv)")
    for f in sorted(glob.glob(os.path.join(OUT, "classify", "*", "*.checks.tsv"))):
        for r in read_tsv(f):
            if r["check"] in ("unjoined", "duplicated_in_output") and r["value"] != "0":
                chk.append({"check": r["check"], "pair": f, "value": r["value"], "expected": 0, "of_cohort": ".",
                            "detail": r["side"]})
                fail.append(f"{f}: {r['check']} {r['value']}")
    write(pd.DataFrame(chk), "checks.tsv")

    # ---- Wilcoxon families ----
    wide = None
    for m in METRICS:
        w = allc.pivot_table(index=["patient", "stratum", "type"], columns="pair", values=m, aggfunc="first")
        w.columns = [f"{m}|{c}" for c in w.columns]
        wide = w if wide is None else wide.join(w)
    wide = wide.reset_index()
    h1 = tests(wide, "H1", [("cross-platform vs cross-method",
                              ["PPR_illumina_aviti", "GATK_illumina_aviti"], ["GATK_PPR_illumina", "GATK_PPR_aviti"])])
    h2 = tests(wide, "H2", [("PPR vs GATK, cross-platform", ["PPR_illumina_aviti"], ["GATK_illumina_aviti"]),
                            ("PPR vs DV, cross-platform", ["PPR_illumina_aviti"], ["DV_illumina_aviti"]),
                            ("DV vs GATK, cross-platform", ["DV_illumina_aviti"], ["GATK_illumina_aviti"])])
    h3 = tests(wide, "H3", [("PPR-DV vs DV-GATK, Illumina", ["DV_PPR_illumina"], ["GATK_DV_illumina"]),
                            ("PPR-DV vs DV-GATK, Aviti", ["DV_PPR_aviti"], ["GATK_DV_aviti"])])
    for name, df in (("wilcoxon_H1.tsv", h1), ("wilcoxon_H2.tsv", h2), ("wilcoxon_H3.tsv", h3)):
        write(df, name)
    print(f"06: {len(g['patient'].unique())} patients, {len(g['pair'].unique())} pairs; H1 {len(h1)} tests, H2 {len(h2)}, "
          f"H3 {len(h3)}; {len(fail)} failed checks")
    for x in fail:
        print("  FAIL " + x)
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
