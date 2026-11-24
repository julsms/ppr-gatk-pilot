#!/usr/bin/env python3
"""05_tables.py - the tables of the unique variants: counts by type and region, carriers and singletons, consequences, novel
alleles, lm in BB, Ti/Tv and het/hom, the GATK status and the RefCall classes (UV_OUT/tables/).
--titv PAIR PATIENT: Ti/Tv and het/hom of one individual and one pair, per side, class and region (UV_OUT/r27/).
Usage: python3 05_tables.py --titv GATK_PPR_illumina <sample>;  python3 05_tables.py
"""
import argparse
import glob
import os
import sys
import time

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "lib"))
import uvlib  # noqa: E402

UV_OUT = os.environ.get("UV_OUT", "uniquevariants")
CONC_OUT = os.environ.get("CONC_OUT", "concordance")
USETS = os.environ.get("USETS", os.path.join(HERE, "config", "unique_sets.tsv"))
PAIRS = os.environ.get("PAIRS", os.path.join(uvlib.CONC_DIR, "config", "pairs.tsv"))
SAMPLE_MAP = os.environ.get("SAMPLE_MAP", os.path.join(uvlib.CONC_DIR, "config", "sample_map.tsv"))
RELATED = os.environ.get("RELATED", "6")
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
CLASS_NAMES = {"GM": "shared", "AM": "gt_discordant", "U": "unique", "NA": "not_assessed"}


def gt_kind(gt):
    al = gt.replace("|", "/").split("/")
    if len(al) == 1:
        return "hom" if al[0] not in ("0", ".") else "other"
    if len(al) == 2 and al[0] == al[1] and al[0] not in ("0", "."):
        return "hom"
    if sorted(al) == ["0", "1"]:
        return "het"
    return "other"


def done_get(path, key):
    if not os.path.exists(path):
        return ""
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 2 and f[0] == key:
            return f[1]
    return ""


def titv(pair, patient, fingerprint=""):
    p = uvlib.tsv_row(PAIRS, "pair", pair)
    c = os.path.join(CONC_OUT, "classify", pair, patient)
    d = os.path.join(UV_OUT, "r27", pair)
    done = os.path.join(d, patient + ".DONE")
    if os.path.exists(done):
        if fingerprint and done_get(done, "fingerprint") == fingerprint:
            print(f"titv {pair} {patient}: already done (current)")
            return
        print(f"STALE\ttitv\t{pair}_{patient}\t{os.path.join(d, patient)}\tthe concordance classification, a prepared file or "
              "the code changed")
        sys.exit(uvlib.EX_STALE)
    if not os.path.exists(c + ".DONE"):
        print(f"MISSING_INPUT\ttitv\t{pair}\t{patient}\t{c}.DONE")
        sys.exit(uvlib.EX_MISSING)
    rec = pd.read_csv(c + ".records.tsv.gz", sep="\t", dtype=str, keep_default_na=False,
                      usecols=["side", "CHROM", "POS", "REF", "ALT", "GT", "TYPE", "class"])
    rec["key"] = uvlib.make_key(rec)
    rows = []
    for s in ("A", "B"):
        st, pl = p[f"{s.lower()}_set"], p[f"{s.lower()}_platform"]
        k = uvlib.load_keys(os.path.join(CONC_OUT, "prep", f"{st}_{pl}", f"{patient}.keys.tsv.gz"))
        strata = uvlib.strata_columns(k)
        r = rec[rec["side"] == s][["key", "class"]].merge(k, on="key", how="inner")
        r["ti"] = [(x, y) in TRANSITIONS for x, y in zip(r["REF"], r["ALT"])] & (r["TYPE"] == "SNV")
        r["tv"] = (r["TYPE"] == "SNV") & ~r["ti"]
        g = r["GT"].map(gt_kind)
        r["het"], r["hom"], r["other_gt"] = g == "het", g == "hom", g == "other"
        r["snv"] = r["TYPE"] == "SNV"
        for cl, name in CLASS_NAMES.items():
            sub = r[r["class"] == cl]
            for stn in ["genome", "autosomes"] + strata:
                m = uvlib.stratum_mask(sub, stn)
                row = {"pair": pair, "patient": patient, "side": s, "set": st, "platform": pl, "class": name, "stratum": stn,
                       "n": int(m.sum())}
                for c2 in ("snv", "ti", "tv", "het", "hom", "other_gt"):
                    row[c2] = int(sub.loc[m, c2].sum())
                rows.append(row)
    os.makedirs(d, exist_ok=True)
    uvlib.write_tsv(pd.DataFrame(rows), os.path.join(d, patient + ".tsv"))
    with open(done + ".tmp", "w") as f:
        f.write(f"date\t{time.strftime('%Y-%m-%d %H:%M:%S')}\nfingerprint\t{fingerprint or '.'}\n"
                f"classify\t{done_get(c + '.DONE', 'fingerprint')}\n")
    os.replace(done + ".tmp", done)


def concat(pattern, done_suffix=None):
    parts = []
    for f in sorted(glob.glob(pattern)):
        if done_suffix and not os.path.exists(f[: -len(done_suffix[0])] + done_suffix[1]):
            continue
        parts.append(pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def num(df, cols):
    for c in cols:
        df[c] = pd.to_numeric(df[c])
    return df


def ratio(a, b):
    return (a / b).where(b > 0).round(6)


def cohort():
    T = os.path.join(UV_OUT, "tables")
    os.makedirs(T, exist_ok=True)
    td = os.path.join(T, "TABLES_DONE")
    if os.path.exists(td):        # the tables are rewritten now: an earlier marker must not vouch for them
        os.replace(td, td + ".previous")
    us = uvlib.read_tsv(USETS)
    patients = [r["patient"] for r in uvlib.read_tsv(SAMPLE_MAP) if r["exclude"] != "yes"]
    checks, fail = [], []

    # unique counts and the checks of extract_unique.py
    uc = concat(os.path.join(UV_OUT, "unique", "*", "*.counts.tsv"), (".counts.tsv", ".DONE"))
    uvlib.write_tsv(uc, os.path.join(T, "unique_counts.tsv"))
    for u in us:
        n = len(glob.glob(os.path.join(UV_OUT, "unique", u["uset"], "*.DONE")))
        checks.append(("unique_people", u["uset"], n))
        for f in glob.glob(os.path.join(UV_OUT, "unique", u["uset"], "*.checks.tsv")):
            if "FAIL" in open(f).read():
                fail.append(f"check of a unique set failed: {f}")

    # DV vs GATK, counts only
    rows = []
    for p in uvlib.read_tsv(PAIRS):
        if p["a_set"] == "GATK" and p["b_set"] == "DV" and p["kind"] == "within":
            for pat in patients:
                c = os.path.join(CONC_OUT, "classify", p["pair"], pat)
                if os.path.exists(c + ".DONE"):
                    rows.append(pd.read_csv(c + ".counts.tsv", sep="\t", dtype=str, keep_default_na=False))
    dv = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if len(dv):
        dv = dv[["patient", "pair", "stratum", "type", "N_A", "N_B", "GM_A", "GM_B", "AM_A", "AM_B", "U_A", "U_B", "LMFLAG_A",
                 "LMFLAG_B", "NA_A", "NA_B"]]
    uvlib.write_tsv(dv, os.path.join(T, "dv_vs_gatk_counts.tsv"))

    # carriers and singletons
    cs = concat(os.path.join(UV_OUT, "carriers", "*.summary.tsv"), (".summary.tsv", ".DONE"))
    uvlib.write_tsv(cs, os.path.join(T, "carrier_summary.tsv"))
    sg = concat(os.path.join(UV_OUT, "carriers", "*.singletons.tsv"), (".singletons.tsv", ".DONE"))
    uvlib.write_tsv(sg, os.path.join(T, "singletons.tsv"))
    if len(sg) and len(cs):
        # primary: summed over the people (each private key is in exactly one person's set), share over the per-person records
        sg = num(sg, ["n", "private_full", "private_full_excl_related", "related"])
        a = sg.groupby(["uset", "stratum", "type"], sort=False)[["n", "private_full"]].sum()
        b = sg[sg["related"] == 0].groupby(["uset", "stratum", "type"], sort=False)[
            ["n", "private_full_excl_related"]].sum().rename(columns={"n": "n_excl_related"})
        ss = a.join(b).reset_index()
        ss["share_private_full"] = ratio(ss["private_full"], ss["n"])
        ss["share_private_full_excl_related"] = ratio(ss["private_full_excl_related"], ss["n_excl_related"])
        # AC = 1 in the multisample VCF, over the distinct keys of the unique set
        c2 = num(cs.copy(), ["n", "in_unique_excl_related", "ac1_unique", "ac1_unique_excl_related"])
        c2 = c2.rename(columns={"n": "n_keys", "in_unique_excl_related": "n_keys_excl_related"})[
            ["uset", "stratum", "type", "n_keys", "ac1_unique", "n_keys_excl_related", "ac1_unique_excl_related"]]
        ss = ss.merge(c2, on=["uset", "stratum", "type"], how="left")
        ss["share_ac1_unique_excl_related"] = ratio(ss["ac1_unique_excl_related"], ss["n_keys_excl_related"])
        ss["share_ac1_unique"] = ratio(ss["ac1_unique"], ss["n_keys"])
        uvlib.write_tsv(ss, os.path.join(T, "singletons_summary.tsv"))
    for f in glob.glob(os.path.join(UV_OUT, "carriers", "*.checks.tsv")):
        if "FAIL" in open(f).read():
            fail.append(f"check of a carrier table failed: {f}")

    # consequences, novel, LM in BB (distinct keys of each unique set, joined to VEP)
    vp = os.path.join(UV_OUT, "vep", "vep_sites.tsv.gz")
    cons, nov, lmbb = [], [], []
    vep = None
    if os.path.exists(vp):
        vep = pd.read_csv(vp, sep="\t", dtype=str, keep_default_na=False,
                          usecols=["key", "consequence_class", "IMPACT", "in_dbSNP", "in_existing_rs", "in_gnomADg",
                                   "in_gnomADe", "in_1000G", "novel", "novel_old"])
        for c in ("in_dbSNP", "in_existing_rs", "in_gnomADg", "in_gnomADe", "in_1000G", "novel", "novel_old"):
            vep[c] = vep[c] == "1"
        for f in open(os.path.join(UV_OUT, "vep", "vep_sites.tsv.gz.checks.tsv")):
            if f.startswith("keys_without_row") and f.split("\t")[1].strip() != "0":
                fail.append("VEP: keys without a row")
    for u in us:
        cf = os.path.join(UV_OUT, "carriers", u["uset"] + ".carriers.tsv.gz")
        if not os.path.exists(cf):
            continue
        k = uvlib.load_keys(cf)
        strata = uvlib.strata_columns(k)
        if vep is not None:
            j = k.merge(vep, on="key", how="left")
            if j["consequence_class"].isna().any():
                fail.append(f"{u['uset']}: {int(j['consequence_class'].isna().sum())} keys without VEP annotation")
            j = j.fillna({"consequence_class": "not_annotated", "IMPACT": "not_annotated"})
            for stn in ["genome", "autosomes"] + strata:
                m = uvlib.stratum_mask(j, stn)
                for (cc, im), n in j[m].groupby(["consequence_class", "IMPACT"]).size().items():
                    cons.append({"uset": u["uset"], "stratum": stn, "consequence_class": cc, "IMPACT": im, "n": int(n),
                                 "n_stratum": int(m.sum())})
            ja = j[j["AUTO"]]
            for stn in ["autosomes"] + strata:
                m = uvlib.stratum_mask(ja, stn)
                sub = ja[m]
                row = {"uset": u["uset"], "stratum": stn, "n": len(sub)}
                for c in ("in_dbSNP", "in_existing_rs", "in_gnomADg", "in_gnomADe", "in_1000G", "novel", "novel_old"):
                    row[c] = int(sub[c].fillna(False).astype(bool).sum())
                row["share_novel"] = round(row["novel"] / row["n"], 6) if row["n"] else float("nan")
                nov.append(row)
        if u["pipeline"] == "PPR" and "LM" in k.columns and "BB" in k.columns:
            n_lm = int(k["LM"].sum())
            lmbb.append({"uset": u["uset"], "patient": "cohort_distinct", "n_LM": n_lm, "n_LM_in_BB": int((k["LM"] & k["BB"]).sum()),
                         "share": round((k["LM"] & k["BB"]).sum() / n_lm, 6) if n_lm else float("nan")})
            for pat in patients:
                kf = os.path.join(UV_OUT, "unique", u["uset"], pat + ".keys.tsv.gz")
                if os.path.exists(kf):
                    # LM is the low-mappability stratum (the lm flag of vcfeval is the column BK_lm)
                    kk = uvlib.load_keys(kf, usecols=["CHROM", "POS", "REF", "ALT", "LM", "BB"])
                    n_lm = int(kk["LM"].sum())
                    lmbb.append({"uset": u["uset"], "patient": pat, "n_LM": n_lm, "n_LM_in_BB": int((kk["LM"] & kk["BB"]).sum()),
                                 "share": round((kk["LM"] & kk["BB"]).sum() / n_lm, 6) if n_lm else float("nan")})
    uvlib.write_tsv(pd.DataFrame(cons), os.path.join(T, "consequences.tsv"))
    uvlib.write_tsv(pd.DataFrame(nov), os.path.join(T, "novel.tsv"))
    uvlib.write_tsv(pd.DataFrame(lmbb), os.path.join(T, "lm_in_bb.tsv"))

    # Ti/Tv and het/hom
    r27 = concat(os.path.join(UV_OUT, "r27", "*", "*.tsv"), (".tsv", ".DONE"))
    if len(r27):
        r27 = num(r27, ["n", "snv", "ti", "tv", "het", "hom", "other_gt"])
        s = r27.groupby(["pair", "platform", "side", "set", "class", "stratum"], sort=False)[
            ["n", "snv", "ti", "tv", "het", "hom", "other_gt"]].sum().reset_index()
        s["people"] = r27.groupby(["pair", "platform", "side", "set", "class", "stratum"], sort=False)["patient"].nunique().values
        s["titv"] = ratio(s["ti"], s["tv"])
        s["hethom"] = ratio(s["het"], s["hom"])
        uvlib.write_tsv(s, os.path.join(T, "r27_titv_hethom.tsv"))
        uvlib.write_tsv(r27, os.path.join(T, "r27_per_person.tsv"))

    # GATK status of the PPR-only alleles, and the RefCall classes
    rc = concat(os.path.join(UV_OUT, "reasons", "*", "*.counts.tsv"), (".counts.tsv", ".DONE"))
    uvlib.write_tsv(rc, os.path.join(T, "r119_classes.tsv"))
    if len(rc):
        cl = [c for c in rc.columns if c[:2] in ("a_", "a2", "a3", "b_", "c_", "d_")]
        rc = num(rc, ["n"] + cl)
        s = rc.groupby(["mode", "platform", "stratum", "type"], sort=False)[["n"] + cl].sum().reset_index()
        for c in cl:
            s["share_" + c] = ratio(s[c], s["n"])
        uvlib.write_tsv(s, os.path.join(T, "r119_summary.tsv"))
    for f in glob.glob(os.path.join(UV_OUT, "reasons", "*", "*.checks.tsv")):
        if "FAIL" in open(f).read():
            fail.append(f"check of the GATK status classes failed: {f}")
    rf = concat(os.path.join(UV_OUT, "refcall", "*", "*.counts.tsv"), (".counts.tsv", ".DONE"))
    uvlib.write_tsv(rf, os.path.join(T, "refcall_classes.tsv"))
    if len(rf):
        rf = num(rf, ["n", "i_refcall", "ii_none"])
        s = rf.groupby(["uset", "stratum", "type"], sort=False)[["n", "i_refcall", "ii_none"]].sum().reset_index()
        s["share_i_refcall"] = ratio(s["i_refcall"], s["n"])
        uvlib.write_tsv(s, os.path.join(T, "refcall_summary.tsv"))

    with open(os.path.join(T, "checks.tsv"), "w") as f:
        f.write("check\ttarget\tvalue\n")
        for c in checks:
            f.write("\t".join(str(x) for x in c) + "\n")
        for x in fail:
            f.write(f"FAIL\t.\t{x}\n")
    if fail:
        sys.exit("05_tables.py: " + "; ".join(fail))
    with open(os.path.join(T, "TABLES_DONE"), "w") as f:
        f.write(time.strftime("%Y-%m-%d %H:%M:%S") + "\n")
    print(f"05: tables written to {T}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--titv", nargs=2, metavar=("PAIR", "PATIENT"))
    ap.add_argument("--fingerprint", default="")
    a = ap.parse_args()
    if a.titv:
        titv(*a.titv, fingerprint=a.fingerprint)
    else:
        cohort()


if __name__ == "__main__":
    main()
