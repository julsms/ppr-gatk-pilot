#!/usr/bin/env python3
"""validation_tables.py - cohort tables of the k-mer support, from the sample_counts.tsv of every sample: counts_long.tsv,
per_sample_genome.tsv, summary.tsv, paired_tests.tsv (paired Wilcoxon signed-rank tests between pipelines, matched-pairs
rank-biserial r, Holm over the regions), hg002.tsv and completeness.tsv; with --vep also the unique records by VEP consequence.
Usage: python3 validation_tables.py --sheet-dir sheet --results results --state sample_state.tsv --out tables [--vep vep_sites.tsv.gz]
"""
import argparse
import hashlib
import math
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from merfin_status import sample_states  # noqa: E402

COMPARISONS = [   # (name, set A, set B): differences are B - A
    ("PPR_ALL-GATK_ALL", "GATK_ALL", "PPR_ALL"),
    ("PPR_ALL-DV_ALL", "DV_ALL", "PPR_ALL"),
    ("DV_ALL-GATK_ALL", "GATK_ALL", "DV_ALL"),
    ("PPRvsGATK-GATKvsPPR", "GATKvsPPR", "PPRvsGATK"),
    ("PPRvsDV-DVvsPPR", "DVvsPPR", "PPRvsDV"),
]


def write(df, path):
    tmp = path + ".tmp"
    df.to_csv(tmp, sep="\t", index=False, na_rep="NA", float_format="%.6g")
    os.replace(tmp, path)


def holm(p):
    """Holm adjustment of a list of p-values (NaN kept)"""
    idx = [i for i, x in enumerate(p) if x == x]
    order = sorted(idx, key=lambda i: p[i])
    m = len(order)
    out = [float("nan")] * len(p)
    run = 0.0
    for k, i in enumerate(order):
        run = max(run, min(1.0, (m - k) * p[i]))
        out[i] = run
    return out


def wilcoxon(d):
    """two-sided Wilcoxon signed-rank on the paired differences d (zeros dropped); p, rank-biserial r"""
    d = [x for x in d if x == x and x != 0]
    if len(d) < 2:
        return float("nan"), float("nan"), len(d)
    try:
        from scipy.stats import rankdata, wilcoxon as wx
    except ImportError:
        return float("nan"), float("nan"), len(d)
    ranks = rankdata([abs(x) for x in d])
    rp = sum(r for r, x in zip(ranks, d) if x > 0)
    rn = sum(r for r, x in zip(ranks, d) if x < 0)
    r = (rp - rn) / (rp + rn)
    p = wx(d).pvalue
    return float(p), float(r), len(d)


def consequence_tables(a, done):
    """the unique records of the done samples by VEP consequence class, with their Merfin status (direct and context)"""
    if not a.vep or not os.path.exists(a.vep):
        print(f"unique by consequence: skipped, no VEP table ({a.vep or 'not given'})")
        return
    vep = pd.read_csv(a.vep, sep="\t", dtype=str, keep_default_na=False, usecols=["key", "consequence_class"])
    vep = vep.set_index("key")["consequence_class"]
    vcfs = pd.read_csv(os.path.join(a.sheet_dir, "vcfs.tsv"), sep="\t", dtype=str, keep_default_na=False)
    rows = []
    for pl, sm, kind in done:
        res = os.path.join(a.results, pl, sm)
        sub = vcfs[(vcfs["platform"] == pl) & (vcfs["sample"] == sm)]
        cache = {}

        def status(sid):
            if sid not in cache:
                cache[sid] = pd.read_csv(os.path.join(res, f"{sid}.status.tsv.gz"), sep="\t", dtype={"key": str, "merfin": int})
            return cache[sid]
        for r in sub[sub["kind"] == "unique"].itertuples(index=False):
            k = pd.read_csv(r.keys, sep="\t", dtype=str, keep_default_na=False, usecols=["CHROM", "POS", "REF", "ALT"])
            keys = set(k["CHROM"] + ":" + k["POS"] + ":" + k["REF"] + ":" + k["ALT"])
            modes = {"context": status(r.context)[status(r.context)["key"].isin(keys)]}
            if r.run == "yes":
                modes["direct"] = status(r.set_id)
            for mode, st in modes.items():
                cons = vep.reindex(st["key"]).fillna("not_annotated").values
                g = st.assign(consequence_class=cons).groupby("consequence_class")["merfin"].agg(["size", "sum"])
                for c, x in g.iterrows():
                    rows.append({"platform": pl, "sample": sm, "sample_kind": kind, "set_id": r.set_id, "mode": mode,
                                 "consequence_class": c, "n": int(x["size"]), "validated": int(x["sum"])})
    per = pd.DataFrame(rows)
    if per.empty:
        return
    per["rate"] = per["validated"] / per["n"]
    write(per, os.path.join(a.out, "unique_consequence_per_sample.tsv"))
    co = per[per["sample_kind"] == "cohort"]
    keys = ["platform", "set_id", "mode", "consequence_class"]
    summ = co.groupby(keys).agg(n_samples=("sample", "nunique"), sum_n=("n", "sum"), sum_validated=("validated", "sum"),
                                median_rate=("rate", "median")).reset_index()
    summ["pooled_rate"] = summ["sum_validated"] / summ["sum_n"]
    write(summ, os.path.join(a.out, "unique_consequence_summary.tsv"))
    print(f"unique by consequence: {len(per)} rows, {per['consequence_class'].nunique()} classes")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet-dir", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--state", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--vep", default="")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    samples = pd.read_csv(os.path.join(a.sheet_dir, "samples.tsv"), sep="\t", dtype=str, keep_default_na=False)
    states = sample_states(samples, a.results, a.state)
    parts, comp, fp, done_list = [], [], hashlib.md5(), []
    for r in samples.itertuples(index=False):
        d = os.path.join(a.results, r.platform, r.sample)
        state = states[(r.platform, r.sample)]
        done = state == "current"
        comp.append({"platform": r.platform, "sample": r.sample, "kind": r.kind,
                     "sample_done": "yes" if os.path.exists(os.path.join(d, "SAMPLE_DONE")) else "no", "state": state})
        if not done:
            continue
        done_list.append((r.platform, r.sample, r.kind))
        fp.update(open(os.path.join(d, "SAMPLE_DONE"), "rb").read())
        c = pd.read_csv(os.path.join(d, "sample_counts.tsv"), sep="\t", dtype={"n": int, "validated": int})
        c.insert(2, "sample_kind", r.kind)
        parts.append(c)
    comp = pd.DataFrame(comp)
    write(comp, os.path.join(a.out, "completeness.tsv"))
    if not parts:
        print("no current sample yet: only completeness.tsv written")
        td = os.path.join(a.out, "TABLES_DONE")
        if os.path.exists(td):
            os.replace(td, td + ".previous")
        return
    long = pd.concat(parts, ignore_index=True)
    long["rate"] = [v / n if n else float("nan") for v, n in zip(long["validated"], long["n"])]
    write(long, os.path.join(a.out, "counts_long.tsv"))
    g = long[(long["stratum"] == "genome") & (long["type"] == "ALL")]
    write(g[["platform", "sample", "sample_kind", "set_id", "kind", "pipeline", "mode", "n", "validated", "rate"]],
          os.path.join(a.out, "per_sample_genome.tsv"))
    write(long[long["sample_kind"] == "hg002"], os.path.join(a.out, "hg002.tsv"))

    co = long[long["sample_kind"] == "cohort"]
    keys = ["platform", "set_id", "kind", "pipeline", "mode", "stratum", "type"]
    rows = []
    for k, sub in co.groupby(keys, sort=False):
        x = sub["rate"].dropna()
        row = dict(zip(keys, k))
        row.update({"n_samples": int(len(x)), "mean": x.mean(), "sd": x.std(), "median": x.median(), "q1": x.quantile(0.25),
                    "q3": x.quantile(0.75), "min": x.min(), "max": x.max(),
                    "pooled_rate": sub["validated"].sum() / sub["n"].sum() if sub["n"].sum() else float("nan"),
                    "median_n": sub["n"].median(), "median_validated": sub["validated"].median()})
        rows.append(row)
    write(pd.DataFrame(rows), os.path.join(a.out, "summary.tsv"))

    tests = []
    piv = co.pivot_table(index=["platform", "mode", "stratum", "type", "sample"], columns="set_id", values="rate")
    for name, sa, sb in COMPARISONS:
        if sa not in piv.columns or sb not in piv.columns:
            continue
        for (pl, mode, typ), sub in piv.reset_index().groupby(["platform", "mode", "type"], sort=False):
            block = []
            for st, s2 in sub.groupby("stratum", sort=False):
                pair = s2[[sa, sb]].dropna()
                if pair.empty:
                    continue
                dif = list(pair[sb] - pair[sa])
                p, r, n_nz = wilcoxon(dif)
                block.append({"platform": pl, "comparison": name, "mode": mode, "stratum": st, "type": typ,
                              "n_pairs": len(pair), "n_nonzero": n_nz, "median_A": pair[sa].median(),
                              "median_B": pair[sb].median(), "median_diff": float(pd.Series(dif).median()),
                              "rank_biserial_r": r, "p": p})
            for row, ph in zip(block, holm([b["p"] for b in block])):
                row["p_holm"] = ph
                if not (math.isnan(row["rank_biserial_r"]) or -1 <= row["rank_biserial_r"] <= 1):
                    sys.exit(f"rank-biserial r out of [-1, 1]: {row}")
            tests += block
    write(pd.DataFrame(tests), os.path.join(a.out, "paired_tests.tsv"))
    consequence_tables(a, done_list)

    n_coh = int((comp["kind"] == "cohort").sum())
    n_done = int(((comp["kind"] == "cohort") & (comp["state"] == "current")).sum())
    other = comp[(comp["state"] != "current") & (comp["sample_done"] == "yes")]
    if len(other):
        print(f"left out (SAMPLE_DONE present, not current): {', '.join(other['platform'] + ':' + other['sample'] + ' ' + other['state'])}")
    print(f"cohort samples current {n_done}/{n_coh}; hg002 current {int(((comp['kind'] == 'hg002') & (comp['state'] == 'current')).sum())}")
    td = os.path.join(a.out, "TABLES_DONE")
    if n_done == n_coh:
        with open(td + ".tmp", "w") as fh:
            fh.write(f"date\t{pd.Timestamp.utcnow():%Y-%m-%d %H:%M:%S} UTC\nfingerprint\t{fp.hexdigest()}\nsamples\t{n_done}\n")
        os.replace(td + ".tmp", td)
        print("TABLES_DONE")
    elif os.path.exists(td):
        os.replace(td, td + ".previous")


if __name__ == "__main__":
    main()
