#!/usr/bin/env python3
"""hg002_merfin_giab.py - Merfin support against the GIAB status of the same HG002 records: merfin_hg002.tsv, merfin_hg002_direct.tsv,
merfin_giab_2x2.tsv and merfin_giab_summary.tsv (records, Merfin share, GIAB precision, their gap, and the share of the GIAB false
positives with k-mer support). Usage: python3 hg002_merfin_giab.py --sheet-dir sheet --results results --state sample_state.tsv
--sets config/merfin_sets.tsv --uv-out <unique-variant output folder> --out hg002
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from merfin_status import sample_states  # noqa: E402


def write(df, path):
    tmp = path + ".tmp"
    df.to_csv(tmp, sep="\t", index=False, na_rep="NA", float_format="%.6g")
    os.replace(tmp, path)


def status(path):
    return pd.read_csv(path, sep="\t", dtype={"key": str, "merfin": int}, keep_default_na=False)


def keys_of(path):
    k = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, usecols=["CHROM", "POS", "REF", "ALT"])
    return set(k["CHROM"] + ":" + k["POS"] + ":" + k["REF"] + ":" + k["ALT"])


def summary_rows(j, ex):
    ev = j[j["giab"].isin(["TP", "FP", "FP_gt"])]
    tp = int((ev["giab"] == "TP").sum())
    fpx = ev[ev["giab"] != "TP"]
    row = dict(ex)
    row.update({"records": len(j), "merfin_share": j["merfin"].mean() if len(j) else float("nan"),
                "evaluable": len(ev), "giab_precision": tp / len(ev) if len(ev) else float("nan"),
                "merfin_share_evaluable": ev["merfin"].mean() if len(ev) else float("nan"),
                "giab_fp": len(fpx), "merfin_supports_giab_fp": fpx["merfin"].mean() if len(fpx) else float("nan")})
    row["gap_merfin_minus_precision"] = row["merfin_share_evaluable"] - row["giab_precision"]
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet-dir", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--state", default="")
    ap.add_argument("--sets", required=True)
    ap.add_argument("--uv-out", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sets = pd.read_csv(a.sets, sep="\t", dtype=str, keep_default_na=False).set_index("set_id")
    samples = pd.read_csv(os.path.join(a.sheet_dir, "samples.tsv"), sep="\t", dtype=str, keep_default_na=False)
    vcfs = pd.read_csv(os.path.join(a.sheet_dir, "vcfs.tsv"), sep="\t", dtype=str, keep_default_na=False)
    hg = samples[samples["kind"] == "hg002"]
    states = sample_states(hg, a.results, a.state)
    full, direct, x2, summ = [], [], [], []
    n_done, truth_missing = 0, []
    for s in hg.itertuples(index=False):
        res = os.path.join(a.results, s.platform, s.sample)
        if states[(s.platform, s.sample)] != "current":
            print(f"HG002 {s.platform}: not used, state {states[(s.platform, s.sample)]}")
            continue
        n_done += 1
        rows = vcfs[(vcfs["platform"] == s.platform) & (vcfs["sample"] == s.sample)]
        st = {r.set_id: status(os.path.join(res, f"{r.set_id}.status.tsv.gz"))
              for r in rows.itertuples(index=False) if os.path.exists(os.path.join(res, f"{r.set_id}.status.tsv.gz"))}
        truth = {}
        for r in rows[rows["kind"] == "all"].itertuples(index=False):
            hname = sets.loc[r.set_id, "hg002_name"]
            f = st[r.set_id].assign(set=hname, platform=s.platform)[["set", "platform", "key", "merfin"]]
            full.append(f)
            tp = os.path.join(a.uv_out, "hg002", "truth", f"{hname}_{s.platform}.status.tsv.gz")
            if os.path.exists(tp):
                truth[r.set_id] = pd.read_csv(tp, sep="\t", dtype=str, keep_default_na=False, usecols=["key", "status"])
            else:
                truth_missing.append(tp)
        for r in rows.itertuples(index=False):
            hname_full = sets.loc[r.context, "hg002_name"] if r.kind == "unique" else sets.loc[r.set_id, "hg002_name"]
            uset = f"{sets.loc[r.set_id, 'hg002_name']}_{s.platform}" if r.kind == "unique" else f"all_{hname_full}"
            base = r.context if r.kind == "unique" else r.set_id
            keys = keys_of(r.keys)
            if r.kind == "unique" and r.set_id in st:
                d = st[r.set_id].assign(uset=uset, set=hname_full, platform=s.platform)[["uset", "set", "platform", "key", "merfin"]]
                direct.append(d)
            if base not in truth:
                continue
            modes = {"context": st[base][st[base]["key"].isin(keys)]}
            if r.kind == "unique" and r.set_id in st:
                modes["direct"] = st[r.set_id]
            for mode, m in modes.items():
                j = m.merge(truth[base], on="key", how="left")
                j["giab"] = j["status"].fillna("NE")
                for (mv, gv), n in j.groupby(["merfin", "giab"]).size().items():
                    x2.append({"uset": uset, "platform": s.platform, "mode": mode, "merfin": int(mv), "giab": gv, "n": int(n)})
                summ.append(summary_rows(j, {"uset": uset, "platform": s.platform, "mode": mode}))
    cols_full = ["set", "platform", "key", "merfin"]
    write(pd.concat(full, ignore_index=True) if full else pd.DataFrame(columns=cols_full), os.path.join(a.out, "merfin_hg002.tsv"))
    write(pd.concat(direct, ignore_index=True) if direct else pd.DataFrame(columns=["uset"] + cols_full),
          os.path.join(a.out, "merfin_hg002_direct.tsv"))
    if x2:
        write(pd.DataFrame(x2), os.path.join(a.out, "merfin_giab_2x2.tsv"))
        write(pd.DataFrame(summ), os.path.join(a.out, "merfin_giab_summary.tsv"))
    else:
        note = pd.DataFrame([{"uset": "truth_not_available", "platform": "NA", "mode": "NA", "merfin": "NA", "giab": "NA", "n": 0}])
        write(note, os.path.join(a.out, "merfin_giab_2x2.tsv"))
    for t in sorted(set(truth_missing)):
        print(f"truth status not available yet (analyses/uniquevariants/hg002_precision.sh): {t}")
    print(f"HG002 samples current {n_done}/{len(hg)}; 2x2 rows {len(x2)}")
    d = os.path.join(a.out, "HG002_MERFIN_DONE")
    if len(hg) and n_done == len(hg):
        with open(d + ".tmp", "w") as fh:
            fh.write(f"samples\t{n_done}\ntruth_missing\t{len(set(truth_missing))}\n")
        os.replace(d + ".tmp", d)
    elif os.path.exists(d):
        os.replace(d, d + ".previous")


if __name__ == "__main__":
    main()
