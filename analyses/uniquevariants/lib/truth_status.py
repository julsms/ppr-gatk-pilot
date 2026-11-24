#!/usr/bin/env python3
"""truth_status.py - the GIAB status of every record of a prepared HG002 call set, and the precision of the unique calls.
status: per record TP, FP, FP_gt (right allele, wrong genotype) inside the benchmark regions (stratum BB of the key tables),
NE outside them; per truth allele TP, GT or NF. tables: precision = TP / (TP + FP + FP_gt) by unique set, region and type, and
check D (D1: shared TP records must be shared in the pair run, at least 0.995; D2 is reported only).
"""
import argparse
import io
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uvlib  # noqa: E402

D1_MIN = 0.995
D2_TOL = 0.02


def query_output(bcftools, vdir):
    res = subprocess.run([bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%BD\t%BK]\n",
                          os.path.join(vdir, "output.vcf.gz")], check=True, stdout=subprocess.PIPE)
    return pd.read_csv(io.StringIO(res.stdout.decode()), sep="\t", header=None, dtype=str, keep_default_na=False, quoting=3,
                       names=["CHROM", "POS", "REF", "ALT", "GT_T", "BD_T", "BK_T", "GT_Q", "BD_Q", "BK_Q"])


def bb_flags(keys, bb, what):
    if bb not in keys.columns:
        sys.exit(f"truth_status.py: no stratum {bb} in the {what} key table (the benchmark regions must be a stratum of the "
                 f"concordance config/strata.tsv)")
    return keys[bb].astype(bool)


def status(a):
    c04 = uvlib.conc04()
    keys = uvlib.load_keys(a.keys, usecols=["CHROM", "POS", "REF", "ALT", "TYPE", a.bb_stratum])
    keys["BB"] = bb_flags(keys, a.bb_stratum, "query")
    tkeys = uvlib.load_keys(a.truth_keys, usecols=["CHROM", "POS", "REF", "ALT", a.bb_stratum])
    tkeys["BB"] = bb_flags(tkeys, a.bb_stratum, "truth")
    known, tknown = set(keys["key"]), set(tkeys["key"])
    out = query_output(a.bcftools, a.vcfeval)
    sq, _, _, unjoined = c04.side_records(out, "Q", known)
    if unjoined:
        sys.exit(f"truth_status.py: {len(unjoined)} query alleles of the output are not in {a.keys}")
    st, _, _, unjoined = c04.side_records(out, "T", tknown)
    if unjoined:
        sys.exit(f"truth_status.py: {len(unjoined)} truth alleles of the output are not in {a.truth_keys} (the truth must be "
                 "prepared as the call sets)")
    # query records: status inside BB from BD and BK, NE outside BB or without an output line
    sq = sq.drop_duplicates("key")
    m = keys[["key", "BB"]].merge(sq[["key", "BD", "BK"]], on="key", how="left")
    m["BD"] = m["BD"].fillna("MISSING")
    m["BK"] = m["BK"].fillna(".")
    m["status"] = ["NE" if not bb else ("TP" if bd == "TP" else ("FP_gt" if bd == "FP" and bk == "am" else
                                                                 ("FP" if bd == "FP" else "NE")))
                   for bb, bd, bk in zip(m["BB"], m["BD"], m["BK"])]
    m["BB"] = m["BB"].astype(int)
    uvlib.write_tsv(m, a.out_prefix + ".status.tsv.gz")
    # truth alleles: every prepared truth record, its BD and BK from the output, its state
    n_dup = int(st.duplicated("key").sum())
    st1 = st.drop_duplicates("key")
    t = tkeys[["key", "BB"]].merge(st1[["key", "BD", "BK"]], on="key", how="left").rename(columns={"key": "tkey"})
    t["BD"] = t["BD"].fillna("MISSING")
    t["BK"] = t["BK"].fillna(".")
    t["state"] = ["TP" if bd == "TP" else ("GT" if bd == "FN" and bk == "am" else "NF") for bd, bk in zip(t["BD"], t["BK"])]
    t["BB"] = t["BB"].astype(int)
    uvlib.write_tsv(t, a.out_prefix + ".truth.tsv.gz")
    # the query records on one output line with a single truth allele, both TP
    single = st[st.groupby("line")["key"].transform("size") == 1]
    qmap = (sq[sq["BD"] == "TP"][["key", "line"]]
            .merge(single[single["BD"] == "TP"][["key", "line"]].rename(columns={"key": "tkey"}), on="line")[["key", "tkey"]])
    uvlib.write_tsv(qmap, a.out_prefix + ".qmap.tsv.gz")
    sc = m["status"].value_counts()
    tc = t.loc[t["BB"] == 1, "state"].value_counts()
    miss_bb = int(((t["BD"] == "MISSING") & (t["BB"] == 1)).sum())
    print("truth status: " + ", ".join(f"{s} {int(sc.get(s, 0))}" for s in ("TP", "FP", "FP_gt", "NE"))
          + f" (query records in BB {int(m['BB'].sum())} of {len(m)}); truth alleles in BB {int(t['BB'].sum())} of {len(t)}: "
          + ", ".join(f"{s} {int(tc.get(s, 0))}" for s in ("TP", "GT", "NF"))
          + f"; truth alleles in BB without an output line {miss_bb}; truth alleles on several lines {n_dup}; qmap {len(qmap)}")


def prec_rows(df, strata, extra, by_cons=False):
    rows = []
    groups = [("ALL", df)]
    if by_cons:
        groups += [(c, df[df["consequence_class"] == c]) for c in sorted(df["consequence_class"].unique())]
    for cons, g in groups:
        for st in (["genome"] if by_cons else ["genome", "autosomes"] + list(strata)):
            m = uvlib.stratum_mask(g, st)
            for t in uvlib.TYPES:
                mm = m if t == "ALL" else m & (g["TYPE"] == t)
                s = g.loc[mm, "status"].value_counts()
                tp, fp, fg, ne = (int(s.get(x, 0)) for x in ("TP", "FP", "FP_gt", "NE"))
                ev = tp + fp + fg
                row = dict(extra)
                row.update({"consequence": cons, "stratum": st, "type": t, "n": int(mm.sum()), "TP": tp, "FP": fp,
                            "FP_gt": fg, "NE": ne,
                            "precision": round(tp / ev, 8) if ev else float("nan"),
                            "FDR": round(1 - tp / ev, 8) if ev else float("nan"),
                            "precision_allele": round((tp + fg) / ev, 8) if ev else float("nan")})
                rows.append(row)
    return rows


def tables(a):
    H = a.hg002_out
    usets = [u for u in uvlib.read_tsv(a.usets) if u["platform"] == a.platform]
    pairs = {p["pair"]: p for p in uvlib.read_tsv(os.path.join(a.pairs))}
    vep = None
    if a.vep and os.path.exists(a.vep):
        vep = pd.read_csv(a.vep, sep="\t", dtype=str, keep_default_na=False, usecols=["key", "consequence_class"])
    merfin = None
    if a.merfin and os.path.exists(a.merfin):
        merfin = pd.read_csv(a.merfin, sep="\t", dtype=str, keep_default_na=False)
    prec, precc, m2 = [], [], []
    done_sets = set()

    def side_table(st, pl):
        keys = uvlib.load_keys(os.path.join(H, "prep", f"{st}_{pl}", "HG002.keys.tsv.gz"))
        s = pd.read_csv(os.path.join(H, "truth", f"{st}_{pl}.status.tsv.gz"), sep="\t", dtype=str, keep_default_na=False,
                        usecols=["key", "status"])
        t = keys.merge(s, on="key", how="left")
        t["status"] = t["status"].fillna("NE")
        t["consequence_class"] = (t[["key"]].merge(vep, on="key", how="left")["consequence_class"].fillna("not_annotated").values
                                  if vep is not None else "not_annotated")
        return t

    for u in usets:
        pr = u["pairs"].split(",")
        s = u["side"].lower()
        st, pl = pairs[pr[0]][f"{s}_set"], pairs[pr[0]][f"{s}_platform"]
        t = side_table(st, pl)
        strata = uvlib.strata_columns(t)
        if (st, pl) not in done_sets:
            done_sets.add((st, pl))
            prec += prec_rows(t, strata, {"uset": f"all_{st}", "platform": pl, "set": st, "class": "all"})
        uk = uvlib.load_keys(os.path.join(H, "unique", u["uset"], "HG002.keys.tsv.gz"), usecols=["CHROM", "POS", "REF", "ALT"])
        tu = t[t["key"].isin(set(uk["key"]))]
        ex = {"uset": u["uset"], "platform": pl, "set": st, "class": "unique"}
        prec += prec_rows(tu, strata, ex)
        precc += prec_rows(tu, strata, ex, by_cons=True)
        if merfin is not None:
            mf = merfin[(merfin["set"] == st) & (merfin["platform"] == pl)][["key", "merfin"]]
            j = tu.merge(mf, on="key", how="left")
            j["merfin"] = j["merfin"].replace("", "NA").fillna("NA")
            for (mv, sv), n in j.groupby(["merfin", "status"]).size().items():
                m2.append({"uset": u["uset"], "platform": pl, "merfin": mv, "giab": sv, "n": int(n)})

    # check D per pair, inside BB (the query status is NE outside BB; the truth alleles are taken inside BB)
    drows = []
    for pname, p in pairs.items():
        if p["a_platform"] != a.platform:
            continue
        rec = pd.read_csv(os.path.join(H, "classify", pname, "HG002.records.tsv.gz"), sep="\t", dtype=str,
                          keep_default_na=False, usecols=["side", "CHROM", "POS", "REF", "ALT", "class"])
        rec["key"] = uvlib.make_key(rec)
        tstate, st_, qm = {}, {}, {}
        for sd in ("a", "b"):
            sset = p[f"{sd}_set"]
            pre = os.path.join(H, "truth", f"{sset}_{a.platform}")
            tr = pd.read_csv(pre + ".truth.tsv.gz", sep="\t", dtype=str, keep_default_na=False, usecols=["tkey", "BB", "state"])
            tr = tr[tr["BB"] == "1"]
            tstate[sd] = dict(zip(tr["tkey"], tr["state"]))
            st_[sd] = pd.read_csv(pre + ".status.tsv.gz", sep="\t", dtype=str, keep_default_na=False, usecols=["key", "status"])
            qm[sd] = pd.read_csv(pre + ".qmap.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
        if set(tstate["a"]) != set(tstate["b"]):
            sys.exit(f"truth_status.py tables: the two truth runs of {pname} do not hold the same truth alleles in BB")
        tpa = {k for k, s in tstate["a"].items() if s == "TP"}
        tpb = {k for k, s in tstate["b"].items() if s == "TP"}
        both = tpa & tpb
        for sd, od, own_tp in (("a", "b", tpa), ("b", "a", tpb)):
            S = sd.upper()
            r = rec[rec["side"] == S]
            shared = set(r.loc[r["class"].isin(["GM", "AM"]), "key"])
            uniq = set(r.loc[r["class"] == "U", "key"])
            q = qm[sd][qm[sd]["tkey"].isin(both)]
            q = q[q["key"].isin(set(st_[sd].loc[st_[sd]["status"] == "TP", "key"]))]
            d1n, d1s = q["key"].nunique(), q.loc[q["key"].isin(shared), "key"].nunique()
            u_tp = int(st_[sd][st_[sd]["key"].isin(uniq) & (st_[sd]["status"] == "TP")].shape[0])
            other = tstate[od]
            only = sum(1 for k in own_tp if other[k] == "NF")
            other_gt = sum(1 for k in own_tp if other[k] == "GT")
            d1 = d1s / d1n if d1n else float("nan")
            d2 = abs(u_tp - only) / only if only else (0.0 if u_tp == 0 else float("inf"))
            drows.append({"pair": pname, "platform": a.platform, "side": S, "set": p[f"{sd}_set"],
                          "truth_TP_this_side": len(own_tp), "truth_TP_both": len(both),
                          "D1_query_TP_on_both": d1n, "D1_shared": d1s,
                          "D1_share": round(d1, 6) if d1 == d1 else float("nan"),
                          "D1_ok": int(d1 != d1 or d1 >= D1_MIN),
                          "D2_unique_TP": u_tp, "D2_truth_only_this_side": only, "D2_truth_other_GT": other_gt,
                          "D2_rel_diff": round(d2, 6), "D2_ok": int(d2 <= D2_TOL),
                          "D2_verdict": "OK" if d2 <= D2_TOL else "WARN"})
    os.makedirs(os.path.join(H, "tables"), exist_ok=True)
    sfx = f"_{a.platform}.tsv"
    uvlib.write_tsv(pd.DataFrame(prec), os.path.join(H, "tables", "precision" + sfx))
    uvlib.write_tsv(pd.DataFrame(precc), os.path.join(H, "tables", "precision_consequence" + sfx))
    uvlib.write_tsv(pd.DataFrame(drows), os.path.join(H, "tables", "check_D" + sfx))
    mf = pd.DataFrame(m2, columns=["uset", "platform", "merfin", "giab", "n"])
    if merfin is None:
        mf = pd.DataFrame([{"uset": "not_available", "platform": a.platform, "merfin": "NA", "giab": "NA", "n": 0}])
    uvlib.write_tsv(mf, os.path.join(H, "tables", "merfin_2x2" + sfx))
    bad1 = [r for r in drows if not r["D1_ok"]]
    warn2 = [r for r in drows if r["D2_verdict"] == "WARN"]
    print(f"tables {a.platform}: {len(prec)} precision rows; check D1: {len(drows) - len(bad1)}/{len(drows)} sides at or above "
          f"{D1_MIN}; check D2: {len(drows) - len(warn2)}/{len(drows)} sides within {D2_TOL}")
    for r in drows:
        print(f"  check D2 {r['D2_verdict']} {r['pair']} side {r['side']} ({r['set']}): unique records GIAB TP in BB "
              f"{r['D2_unique_TP']}, truth alleles in BB found only by this side {r['D2_truth_only_this_side']} "
              f"(other side with another genotype {r['D2_truth_other_GT']}), relative difference {r['D2_rel_diff']}")
    if warn2:
        print(f"WARNING: check D2 above {D2_TOL} for " + ", ".join(f"{r['pair']} side {r['side']}" for r in warn2)
              + f" (reported, not a gate; see check_D_{a.platform}.tsv)")
    if bad1:
        sys.exit("truth_status.py tables: check D1 failed for " + ", ".join(f"{r['pair']} side {r['side']}" for r in bad1)
                 + f" (D1 share >= {D1_MIN}; see check_D_{a.platform}.tsv)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["status", "tables"])
    ap.add_argument("--vcfeval")
    ap.add_argument("--keys")
    ap.add_argument("--truth-keys")
    ap.add_argument("--bb-stratum", default=os.environ.get("BB_STRATUM", "BB"))
    ap.add_argument("--out-prefix")
    ap.add_argument("--usets")
    ap.add_argument("--pairs")
    ap.add_argument("--platform")
    ap.add_argument("--hg002-out")
    ap.add_argument("--vep")
    ap.add_argument("--merfin")
    ap.add_argument("--bcftools", default=os.environ.get("BCFTOOLS", "bcftools"))
    a = ap.parse_args()
    status(a) if a.mode == "status" else tables(a)


if __name__ == "__main__":
    main()
