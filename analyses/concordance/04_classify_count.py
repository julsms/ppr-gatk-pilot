#!/usr/bin/env python3
"""Classifies every record of both call sets of one vcfeval run (GA4GH output; A = TRUTH column, B = QUERY column) and counts the
classes per region and variant type: GM same allele and genotype (BD=TP), AM same allele with another genotype (BD=FP or FN, BK=am),
U unique (BD=FP or FN, BK=lm or none), NA not assessed. J_al = Sbar / (Sbar + U_A + U_B), J_gt = Gbar / (Sbar + U_A + U_B).
Usage: python3 04_classify_count.py --pair <pair> --patient <individual>    (in the output folder, or with OUT set)
"""
import argparse
import csv
import datetime
import hashlib
import io
import os
import resource
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "lib"))
import strata_names  # noqa: E402  (lib/strata_names.py: the stratum name rule and the column prefix)

TYPES = ["ALL", "SNV", "INDEL", "MNP", "OTHER"]
SIDE_COLS = ["N", "GM", "AM", "U", "LMFLAG", "NA", "NA_N", "NA_MISSING", "SPANDEL", "SPANDEL_GM", "SPANDEL_AM", "SPANDEL_U",
             "SPANDEL_NA", "SPANDEL_PAIR", "SPANDEL_PAIR_GM", "SPANDEL_PAIR_AM", "SPANDEL_PAIR_U", "SPANDEL_PAIR_NA"]
NON_STRATA = ("CHROM", "POS", "REF", "ALT", "GT", "TYPE", "SPANDEL", "AUTO", "key")
NA_ALARM = 0.001
EX_MISSING, EX_STALE = 64, 65


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


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


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def config_paths(pair, patient):
    """the paths of one pair and patient under OUT"""
    out = os.environ.get("OUT", ".")
    pairs = os.environ.get("PAIRS", os.path.join(HERE, "config", "pairs.tsv"))
    row = [r for r in read_tsv(pairs) if r["pair"] == pair]
    if not row:
        sys.exit(f"04_classify_count.py: unknown pair {pair}")
    r = row[0]
    a_dir = os.path.join(out, "prep", f"{r['a_set']}_{r['a_platform']}")
    b_dir = os.path.join(out, "prep", f"{r['b_set']}_{r['b_platform']}")
    a = os.path.join(a_dir, f"{patient}.keys.tsv.gz")
    b = os.path.join(b_dir, f"{patient}.keys.tsv.gz")
    v = os.path.join(out, "vcfeval", pair, patient)
    o = os.path.join(out, "classify", pair, patient)
    up = {"vcfeval_fp": done_get(v + ".DONE", "fingerprint"),
          "a_prep_fp": done_get(os.path.join(a_dir, f"{patient}.DONE"), "fingerprint"),
          "b_prep_fp": done_get(os.path.join(b_dir, f"{patient}.DONE"), "fingerprint")}
    return a, b, v, o, r, up


def trim(pos, ref, alt):
    """minimal allele pair: common suffix, then common prefix, keeping one base on each side"""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return pos, ref, alt


def load_keys(path):
    """a key table of 01 with the 0/1 columns as booleans and every stratum column renamed stratum:<name>; and the strata"""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    strata = strata_names.check([c for c in df.columns if c not in NON_STRATA], f"04_classify_count.py: {path}")
    df["key"] = df["CHROM"] + ":" + df["POS"] + ":" + df["REF"] + ":" + df["ALT"]
    for c in df.columns:
        if c not in ("CHROM", "POS", "REF", "ALT", "GT", "TYPE", "key"):
            df[c] = df[c] == "1"        # SPANDEL, AUTO and the strata flags
    return strata_names.namespaced(df, strata), strata


def side_records(out_df, col, known):
    """the (key, BD, BK, output line, BS) of one side of the vcfeval output, joined to the prepared keys"""
    rows, n_exact, n_trim, unjoined = [], 0, 0, []
    gts, bds, bks = out_df[f"GT_{col}"], out_df[f"BD_{col}"], out_df[f"BK_{col}"]
    bss = out_df["BS"] if "BS" in out_df.columns else ["."] * len(out_df)     # callers that query no BS still work
    for line, (chrom, pos, ref, alts, gt, bd, bk, bs) in enumerate(zip(out_df["CHROM"], out_df["POS"], out_df["REF"],
                                                                        out_df["ALT"], gts, bds, bks, bss)):
        idx = {int(x) for x in gt.replace("|", "/").split("/") if x not in (".", "") and x != "0"}
        if not idx:
            continue
        alt_list = alts.split(",")
        for i in sorted(idx):
            alt = alt_list[i - 1]
            k = f"{chrom}:{pos}:{ref}:{alt}"
            if k in known:
                n_exact += 1
            else:
                p2, r2, a2 = trim(int(pos), ref, alt)
                k = f"{chrom}:{p2}:{r2}:{a2}"
                if k in known:
                    n_trim += 1
                else:
                    unjoined.append(k)
                    continue
            rows.append((k, bd, bk, line, bs))
    df = pd.DataFrame(rows, columns=["key", "BD", "BK", "line", "BS"])
    return df, n_exact, n_trim, unjoined


def classify(bd, bk):
    if bd == "TP":
        return "GM"
    if bd in ("FP", "FN"):
        return "AM" if bk == "am" else "U"
    return "NA"


def spandel_pairs(sr, sides):
    """the keys of each side flagged SPANDEL_PAIR, and the number added through the superlocus rule"""
    own = {s: set(sides[s].loc[sides[s]["SPANDEL"], "key"]) for s in ("A", "B")}
    lines = {s: set(sr[s]["line"]) for s in ("A", "B")}
    flag_lines = set()
    for s in ("A", "B"):
        flag_lines |= set(sr[s].loc[sr[s]["key"].isin(own[s]), "line"])
    pair = {s: own[s] | set(sr[s].loc[sr[s]["line"].isin(flag_lines), "key"]) for s in ("A", "B")}
    by_bs = {"A": 0, "B": 0}
    for s, o in (("A", "B"), ("B", "A")):
        x = sr[s]
        shared = x["key"].isin(own[s]) & x["cls"].isin(["GM", "AM"]) & ~x["line"].isin(lines[o])
        bs_set = {b for v in x.loc[shared, "BS"] for b in str(v).split(",") if b not in (".", "")}
        if not bs_set:
            continue
        y = sr[o]
        # shared records of the other side in the same superlocus that have no partner of their own on their line
        hit = (y["cls"].isin(["GM", "AM"]) & ~y["line"].isin(lines[s])
               & y["BS"].map(lambda v: any(b in bs_set for b in str(v).split(","))))
        new = set(y.loc[hit, "key"]) - pair[o]
        by_bs[o] += len(new)
        pair[o] |= new
    return pair, by_bs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair")
    ap.add_argument("--patient")
    ap.add_argument("--a-keys")
    ap.add_argument("--b-keys")
    ap.add_argument("--vcfeval")
    ap.add_argument("--out-prefix")
    ap.add_argument("--fingerprint", default="")
    ap.add_argument("--bcftools", default=os.environ.get("BCFTOOLS", "bcftools"))
    a = ap.parse_args()

    meta = {"a_set": ".", "a_platform": ".", "b_set": ".", "b_platform": ".", "kind": "."}
    up = {"vcfeval_fp": ".", "a_prep_fp": ".", "b_prep_fp": "."}
    if a.pair and a.patient and not a.a_keys:
        a_keys, b_keys, vdir, prefix, meta, up = config_paths(a.pair, a.patient)
    else:
        a_keys, b_keys, vdir, prefix = a.a_keys, a.b_keys, a.vcfeval, a.out_prefix
        if a.pair and a.patient:
            meta = config_paths(a.pair, a.patient)[4]
    pair, patient = a.pair or ".", a.patient or "."
    if os.path.exists(prefix + ".DONE"):
        if a.fingerprint and done_get(prefix + ".DONE", "fingerprint") == a.fingerprint:
            print(f"04 {pair} {patient}: already done (current)")
            sys.exit(0)
        print(f"STALE\tclassify\t{pair}_{patient}\t{prefix}\tthe vcfeval output, a prepared file or the code changed")
        sys.exit(EX_STALE)
    if not os.path.exists(os.path.join(vdir, "output.vcf.gz")) or not os.path.exists(vdir + ".DONE"):
        print(f"MISSING_INPUT\tclassify\t{pair}\t{patient}\t{vdir}")
        sys.exit(EX_MISSING)
    os.makedirs(os.path.dirname(prefix), exist_ok=True)

    (ka, st_a), (kb, st_b) = load_keys(a_keys), load_keys(b_keys)
    # the strata of side A are counted (a subset key table of the unique-variant analysis has none); side B must have them all
    if [s for s in st_a if s not in st_b]:
        sys.exit(f"04_classify_count.py: strata of side A missing from side B: {[s for s in st_a if s not in st_b]}")
    sides = {"A": ka, "B": kb}
    strata = ["genome", "autosomes"] + st_a
    res = subprocess.run([a.bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/BS[\t%GT\t%BD\t%BK]\n",
                          os.path.join(vdir, "output.vcf.gz")], check=True, stdout=subprocess.PIPE)
    out_df = pd.read_csv(io.StringIO(res.stdout.decode()), sep="\t", header=None, dtype=str, keep_default_na=False,
                         quoting=3, names=["CHROM", "POS", "REF", "ALT", "BS", "GT_T", "BD_T", "BK_T", "GT_Q", "BD_Q", "BK_Q"])

    checks, fail, sr = [], [], {}
    for s, col in (("A", "T"), ("B", "Q")):
        known = set(sides[s]["key"])
        x, n_exact, n_trim, unjoined = side_records(out_df, col, known)
        dup = int(x["key"].duplicated().sum())
        checks += [(s, "joined_exact", n_exact), (s, "joined_after_trim", n_trim), (s, "unjoined", len(unjoined)),
                   (s, "duplicated_in_output", dup)]
        if unjoined:
            fail.append(f"side {s}: {len(unjoined)} output alleles not found in the prepared file")
        if dup:
            fail.append(f"side {s}: {dup} prepared records appear twice in the output")
            x = x.drop_duplicates("key")
        x["cls"] = [classify(bd, bk) for bd, bk in zip(x["BD"], x["BK"])]
        sr[s] = x
    pair_keys, by_bs = spandel_pairs(sr, sides)

    records = []
    for s in ("A", "B"):
        m = sides[s].merge(sr[s][["key", "BD", "BK"]], on="key", how="left")
        m["BD"] = m["BD"].fillna("MISSING")
        m["BK"] = m["BK"].fillna(".")
        m["class"] = [classify(bd, bk) for bd, bk in zip(m["BD"], m["BK"])]
        m["LMFLAG"] = ((m["class"] == "U") & (m["BK"] == "lm")).astype(int)     # not "LM": that is a stratum name
        m["NA_MISSING"] = (m["BD"] == "MISSING").astype(int)
        m["SPANDEL_PAIR"] = m["key"].isin(pair_keys[s])
        for v, c in sorted(m["BD"].value_counts().to_dict().items()):
            checks.append((s, f"BD_{v}", int(c)))
        checks += [(s, "spandel_own", int(m["SPANDEL"].sum())), (s, "spandel_pair", int(m["SPANDEL_PAIR"].sum())),
                   (s, "spandel_pair_by_superlocus", by_bs[s])]
        sides[s] = m
        rec = m[["CHROM", "POS", "REF", "ALT", "GT", "TYPE", "SPANDEL", "SPANDEL_PAIR", "class", "BD", "BK"]].copy()
        rec["SPANDEL"] = rec["SPANDEL"].astype(int)
        rec["SPANDEL_PAIR"] = rec["SPANDEL_PAIR"].astype(int)
        rec.insert(0, "side", s)
        records.append(rec)

    # counts per stratum x type x side
    rows = []
    for st in strata:
        for t in TYPES:
            row = {"patient": patient, "pair": pair, **{x: meta.get(x, ".") for x in ("a_set", "a_platform", "b_set",
                                                                                     "b_platform", "kind")},
                   "stratum": st, "type": t}
            for s in ("A", "B"):
                m = sides[s]
                mask = pd.Series(True, index=m.index)
                if st == "autosomes":
                    mask &= m["AUTO"]
                elif st != "genome":
                    mask &= m[strata_names.COL_PREFIX + st]
                if t != "ALL":
                    mask &= m["TYPE"] == t
                sub = m[mask]
                cls = sub["class"].value_counts()
                spc = sub.loc[sub["SPANDEL"], "class"].value_counts()
                ppc = sub.loc[sub["SPANDEL_PAIR"], "class"].value_counts()
                v = {"N": len(sub), "GM": int(cls.get("GM", 0)), "AM": int(cls.get("AM", 0)), "U": int(cls.get("U", 0)),
                     "LMFLAG": int(sub["LMFLAG"].sum()), "NA": int(cls.get("NA", 0)),
                     "NA_MISSING": int(sub["NA_MISSING"].sum()), "SPANDEL": int(sub["SPANDEL"].sum()),
                     "SPANDEL_PAIR": int(sub["SPANDEL_PAIR"].sum())}
                for c in ("GM", "AM", "U", "NA"):
                    v[f"SPANDEL_{c}"] = int(spc.get(c, 0))
                    v[f"SPANDEL_PAIR_{c}"] = int(ppc.get(c, 0))
                v["NA_N"] = v["NA"] - v["NA_MISSING"]
                if v["N"] != v["GM"] + v["AM"] + v["U"] + v["NA"]:
                    fail.append(f"closure fails: side {s}, {st}, {t}")
                for c in SIDE_COLS:
                    row[f"{c}_{s}"] = v[c]
            rows.append(row)
    cnt = pd.DataFrame(rows)

    def ratio(x, y):
        return (x / y).where(y > 0)

    for s in ("A", "B"):
        cnt[f"S_{s}"] = cnt[f"GM_{s}"] + cnt[f"AM_{s}"]
    sbar = (cnt["S_A"] + cnt["S_B"]) / 2
    gbar = (cnt["GM_A"] + cnt["GM_B"]) / 2
    den = sbar + cnt["U_A"] + cnt["U_B"]
    cnt["Sbar"], cnt["Gbar"] = sbar, gbar
    cnt["J_al"] = ratio(sbar, den)
    cnt["J_gt"] = ratio(gbar, den)
    cnt["share_A"] = ratio(cnt["S_A"], cnt["N_A"])
    cnt["share_B"] = ratio(cnt["S_B"], cnt["N_B"])
    cnt["gt_agree"] = ratio(gbar, sbar)
    cnt["rep_gap"] = (cnt["S_A"] - cnt["S_B"]).abs()
    cnt["NA_share"] = ratio(cnt["NA_A"] + cnt["NA_B"], cnt["N_A"] + cnt["N_B"])
    cnt["NA_alarm"] = (cnt["NA_share"] > NA_ALARM).astype(int)
    # without SPANDEL_PAIR records on both sides (the flagged record and its partner)
    s_ns = sum(cnt[f"S_{s}"] - cnt[f"SPANDEL_PAIR_GM_{s}"] - cnt[f"SPANDEL_PAIR_AM_{s}"] for s in ("A", "B")) / 2
    g_ns = sum(cnt[f"GM_{s}"] - cnt[f"SPANDEL_PAIR_GM_{s}"] for s in ("A", "B")) / 2
    den_ns = s_ns + sum(cnt[f"U_{s}"] - cnt[f"SPANDEL_PAIR_U_{s}"] for s in ("A", "B"))
    cnt["J_gt_nospandel"] = ratio(g_ns, den_ns)
    # back-calculation: U_A + U_B from J_al and Sbar must equal the counted U
    back = (sbar / cnt["J_al"] - sbar).where(cnt["J_al"] > 0)
    bad_back = ((back - (cnt["U_A"] + cnt["U_B"])).abs() > 1e-6 * (den + 1)) & back.notna()
    if bad_back.any():
        fail.append(f"back-calculation of U from J_al fails in {int(bad_back.sum())} rows")
    for c in ["J_al", "J_gt", "share_A", "share_B", "gt_agree", "NA_share", "J_gt_nospandel", "Sbar", "Gbar"]:
        cnt[c] = cnt[c].round(8)

    g = cnt[(cnt["stratum"] == "genome") & (cnt["type"] == "ALL")].iloc[0]
    checks += [(".", "closure", "ok" if not [f for f in fail if f.startswith("closure")] else "FAIL"),
               (".", "NA_share_genome_all", g["NA_share"]), (".", "NA_alarm_genome_all", int(g["NA_alarm"]))]

    cnt.to_csv(prefix + ".counts.tsv.tmp", sep="\t", index=False, na_rep="NA")
    pd.concat(records).to_csv(prefix + ".records.tsv.gz.tmp", sep="\t", index=False,
                              compression={"method": "gzip", "mtime": 0})
    with open(prefix + ".checks.tsv.tmp", "w") as f:
        f.write("side\tcheck\tvalue\n")
        for c in checks:
            f.write("\t".join(str(x) for x in c) + "\n")
        for x in fail:
            f.write(f".\tFAIL\t{x}\n")
    if fail:
        os.replace(prefix + ".checks.tsv.tmp", prefix + ".checks.tsv")
        sys.exit("04_classify_count.py: " + "; ".join(fail))
    for ext in (".counts.tsv", ".records.tsv.gz", ".checks.tsv"):
        os.replace(prefix + ext + ".tmp", prefix + ext)
    rss_mb = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss + resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss) / 1024
    done = [("date", datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")),
            ("fingerprint", a.fingerprint or "."), *up.items(),
            ("counts_md5", md5(prefix + ".counts.tsv")), ("records_md5", md5(prefix + ".records.tsv.gz")),
            ("peak_rss_mb", f"{rss_mb:.0f}"),
            ("genome_all", f"N_A {int(g['N_A'])} N_B {int(g['N_B'])} J_al {g['J_al']} J_gt {g['J_gt']} NA_share {g['NA_share']}")]
    with open(prefix + ".DONE.tmp", "w") as f:
        f.write("".join(f"{k}\t{v}\n" for k, v in done))
    os.replace(prefix + ".DONE.tmp", prefix + ".DONE")
    print(f"04 {pair} {patient}: genome ALL N_A {int(g['N_A'])} N_B {int(g['N_B'])} GM {int(g['GM_A'])}/{int(g['GM_B'])} "
          f"AM {int(g['AM_A'])}/{int(g['AM_B'])} U {int(g['U_A'])}/{int(g['U_B'])} NA {int(g['NA_A'])}/{int(g['NA_B'])} "
          f"J_al {g['J_al']:.4f} J_gt {g['J_gt']:.4f}; peak RSS {rss_mb:.0f} MB")


if __name__ == "__main__":
    main()
