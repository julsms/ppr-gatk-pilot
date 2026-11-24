#!/usr/bin/env python3
"""02_carrier_table.py - carriers, singletons and the multisample VCF of one unique set, over the cohort.
Per distinct allele: the individuals that carry it in the full call set (n_carriers_full), in how many it is unique, shared or
not assessed, the singleton flag (private_full: one carrier) with and without RELATED, and the allele copies (ac_unique).
Usage: python3 02_carrier_table.py --uset PPRvsGATK_illumina
"""
import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import time

import numpy as np
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
BCFTOOLS = os.environ.get("BCFTOOLS", "bcftools")
CODE = ("02_carrier_table.py", "lib/uvlib.py")   # the code that shapes the outputs, relative to this folder


def allele_copies(gt):
    """copies of the ALT allele of a split record: 0/1, ./1 and a haploid 1 count 1, 1/1 counts 2"""
    return sum(1 for x in gt.replace("|", "/").split("/") if x not in ("0", ".", ""))


def code_md5(files):
    """md5 of "name md5" lines of the code files"""
    lines = []
    for f in files:
        with open(os.path.join(HERE, f), "rb") as fh:
            lines.append(f"{f} {hashlib.md5(fh.read()).hexdigest()}\n")
    return hashlib.md5("".join(lines).encode()).hexdigest()


def signature(files):
    h = hashlib.md5()
    for f in files:
        h.update(f.encode())
        h.update(open(f).read().encode() if os.path.exists(f) else b"MISSING")
    h.update(f"related {RELATED}".encode())
    h.update(f"code {code_md5(CODE)}".encode())
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--uset", required=True)
    a = ap.parse_args()
    us = uvlib.tsv_row(USETS, "uset", a.uset)
    pairs = us["pairs"].split(",")
    s = us["side"].lower()
    prow = [uvlib.tsv_row(PAIRS, "pair", p) for p in pairs]
    pset, pplat = prow[0][f"{s}_set"], prow[0][f"{s}_platform"]
    patients = [r["patient"] for r in uvlib.read_tsv(SAMPLE_MAP) if r["exclude"] != "yes"]

    prep = {p: os.path.join(CONC_OUT, "prep", f"{pset}_{pplat}", p) for p in patients}
    cls = {(pr, p): os.path.join(CONC_OUT, "classify", pr, p) for pr in pairs for p in patients}
    uq = {p: os.path.join(UV_OUT, "unique", a.uset, p) for p in patients}
    full_pat = [p for p in patients if os.path.exists(prep[p] + ".DONE")]
    uset_pat = [p for p in patients if os.path.exists(uq[p] + ".DONE")]

    D = os.path.join(UV_OUT, "carriers")
    os.makedirs(D, exist_ok=True)
    P = os.path.join(D, a.uset)
    sig = signature([prep[p] + ".DONE" for p in patients] + [c + ".DONE" for c in cls.values()]
                    + [uq[p] + ".DONE" for p in patients])
    if os.path.exists(P + ".DONE") and os.path.exists(P + ".signature") and open(P + ".signature").read().strip() == sig:
        print(f"02 {a.uset}: already done for these inputs")
        return
    old = [P + x for x in (".carriers.tsv.gz", ".singletons.tsv", ".summary.tsv", ".multisample.vcf.gz",
                           ".multisample.vcf.gz.tbi", ".checks.tsv", ".signature", ".DONE") if os.path.exists(P + x)]
    if old:
        sd = os.path.join(D, f"_superseded_{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}")
        os.makedirs(sd, exist_ok=True)
        for f in old:
            shutil.move(f, sd)
        print(f"02 {a.uset}: the inputs or the code changed; old outputs moved to {sd}")
    if not uset_pat:
        print(f"MISSING_INPUT\tcarriers\t{a.uset}\tno person has an output of this unique set")
        sys.exit(uvlib.EX_MISSING)

    fail, checks = [], []
    # pass 1: the union of the unique keys, with the attributes of the first person that has the key
    base, pkeys, pcopies, n_unique = [], {}, {}, {}
    known = set()
    for p in uset_pat:
        k = uvlib.load_keys(uq[p] + ".keys.tsv.gz")
        n_unique[p] = len(k)
        new = k[~k["key"].isin(known)]
        base.append(new.drop(columns=["GT", "BK_lm"], errors="ignore"))
        known.update(new["key"])
        pkeys[p] = k["key"]
        pcopies[p] = np.array([allele_copies(g) for g in k["GT"]], dtype=np.int32)
    U = pd.concat(base, ignore_index=True)
    idx = pd.Index(U["key"])
    if not idx.is_unique:
        sys.exit("02: duplicated keys in the union")
    nk = len(U)
    strata = uvlib.strata_columns(U)

    # pass 2: carriers in the full prepared sets
    car = np.zeros(nk, dtype=np.int32)
    car_x = np.zeros(nk, dtype=np.int32)
    for p in full_pat:
        k = pd.read_csv(prep[p] + ".keys.tsv.gz", sep="\t", dtype=str, keep_default_na=False,
                        usecols=["CHROM", "POS", "REF", "ALT"])
        c = idx.get_indexer(uvlib.make_key(k))
        c = c[c >= 0]
        car[c] += 1
        if p != RELATED:
            car_x[c] += 1

    # pass 3: unique / shared / not assessed per person, from the classes of the pair runs
    nU = np.zeros(nk, dtype=np.int32)
    nU_x = np.zeros(nk, dtype=np.int32)
    nS = np.zeros(nk, dtype=np.int32)
    nNA = np.zeros(nk, dtype=np.int32)
    side = us["side"]
    for p in uset_pat:
        st = None
        for pr in pairs:
            r = pd.read_csv(cls[(pr, p)] + ".records.tsv.gz", sep="\t", dtype=str, keep_default_na=False,
                            usecols=["side", "CHROM", "POS", "REF", "ALT", "class"])
            r = r[r["side"] == side]
            r = pd.DataFrame({"key": uvlib.make_key(r).values, "class": r["class"].values}).set_index("key")["class"]
            st = r.to_frame(pr) if st is None else st.join(r.rename(pr), how="outer")
        isU = (st == "U").all(axis=1)
        isS = st.isin(["GM", "AM"]).any(axis=1)
        status = np.where(isU, "U", np.where(isS, "S", "NA"))
        c = idx.get_indexer(st.index)
        m = c >= 0
        c, status = c[m], status[m]
        nU[c[status == "U"]] += 1
        nS[c[status == "S"]] += 1
        nNA[c[status == "NA"]] += 1
        if p != RELATED:
            nU_x[c[status == "U"]] += 1
        if int((status == "U").sum()) != n_unique[p]:
            fail.append(f"{p}: {int((status == 'U').sum())} U records in the pair runs, {n_unique[p]} in the unique set")
    if (car[nU > 0] < nU[nU > 0]).any():
        fail.append("a key is unique in more people than carry it in the full set")
    # allele copies in the unique sets, with and without RELATED
    ac = np.zeros(nk, dtype=np.int32)
    ac_x = np.zeros(nk, dtype=np.int32)
    for p in uset_pat:
        c = idx.get_indexer(pkeys[p])
        np.add.at(ac, c, pcopies[p])
        if p != RELATED:
            np.add.at(ac_x, c, pcopies[p])

    n_out = nU + nS + nNA
    U["n_full_patients"] = len(full_pat)
    U["n_uset_patients"] = len(uset_pat)
    U["n_carriers_full"] = car
    U["n_carriers_full_excl_related"] = car_x
    U["n_U"], U["n_S"], U["n_NA"] = nU, nS, nNA
    U["category"] = np.where(nU == n_out, "unique_all", np.where(nS > 0, "unique_some_shared_some", "unique_some_na_some"))
    U["private_full"] = car == 1
    U["private_full_excl_related"] = car_x == 1
    U["ac_unique"], U["ac_unique_excl_related"] = ac, ac_x
    U["in_unique_excl_related"] = nU_x > 0
    U["ac1_unique"] = ac == 1
    U["ac1_unique_excl_related"] = ac_x == 1
    for c in ("unique_all", "unique_some_shared_some", "unique_some_na_some"):
        U[c] = U["category"] == c

    # per-person singleton table
    flags = ["private_full", "private_full_excl_related", "ac1_unique", "ac1_unique_excl_related"]
    rows = []
    for p in uset_pat:
        sub = U.iloc[idx.get_indexer(pkeys[p])].copy()
        # a flag counts for the person who carries the key (each key of the person's set is the person's own)
        if p == RELATED:
            sub["private_full_excl_related"] = False
            sub["ac1_unique_excl_related"] = False
        t = uvlib.count_by(sub, strata, extra={"uset": a.uset, "patient": p, "related": int(p == RELATED)}, value_cols=flags)
        rows.append(t)
    sing = pd.concat(rows, ignore_index=True)
    summ = uvlib.count_by(U, strata, extra={"uset": a.uset, "n_full_patients": len(full_pat),
                                            "n_uset_patients": len(uset_pat)},
                          value_cols=["unique_all", "unique_some_shared_some", "unique_some_na_some", "in_unique_excl_related"]
                          + flags)

    # multisample VCF
    vcfs = [uq[p] + ".vcf.gz" for p in uset_pat]
    ms = P + ".multisample.vcf.gz"
    if len(vcfs) == 1:
        shutil.copy(vcfs[0], ms + ".tmp")
    else:
        subprocess.run([BCFTOOLS, "merge", "-m", "none", "-Oz", "-o", ms + ".tmp"] + vcfs, check=True)
    n_rec, n_gt, n_copies = 0, 0, 0
    for line in uvlib.bcftools_lines(BCFTOOLS, ["query", "-f", "[%GT\t]\n", ms + ".tmp"]):
        n_rec += 1
        gts = line.rstrip("\t\n").split("\t")
        n_gt += sum(1 for g in gts if any(ch in g for ch in "123456789"))
        n_copies += sum(allele_copies(g) for g in gts)
    os.replace(ms + ".tmp", ms)
    subprocess.run(["tabix", "-f", "-p", "vcf", ms], check=True)
    if n_rec != nk:
        fail.append(f"multisample VCF: {n_rec} records, {nk} distinct keys")
    if n_gt != sum(n_unique.values()):
        fail.append(f"multisample VCF: {n_gt} carried genotypes, {sum(n_unique.values())} per-person unique records")
    if int(nU.sum()) != sum(n_unique.values()):
        fail.append(f"n_U adds up to {int(nU.sum())}, the per-person sets to {sum(n_unique.values())}")
    if n_copies != int(ac.sum()):
        fail.append(f"multisample VCF: {n_copies} allele copies, ac_unique adds up to {int(ac.sum())}")

    checks += [("full_patients", len(full_pat)), ("uset_patients", len(uset_pat)),
               ("missing_uset_patients", ",".join(p for p in patients if p not in uset_pat) or "none"),
               ("distinct_keys", nk), ("sum_unique_per_person", sum(n_unique.values())),
               ("multisample_records", n_rec), ("multisample_carried_gt", n_gt), ("sum_n_U", int(nU.sum())),
               ("multisample_allele_copies", n_copies), ("sum_ac_unique", int(ac.sum()))]
    out = U.drop(columns=["key", "unique_all", "unique_some_shared_some", "unique_some_na_some"])
    for c in out.columns:
        if out[c].dtype == bool:
            out[c] = out[c].astype(int)
    uvlib.write_tsv(out, P + ".carriers.tsv.gz")
    uvlib.write_tsv(sing, P + ".singletons.tsv")
    uvlib.write_tsv(summ, P + ".summary.tsv")
    with open(P + ".checks.tsv", "w") as f:
        f.write("check\tvalue\n")
        for c in checks:
            f.write(f"{c[0]}\t{c[1]}\n")
        for x in fail:
            f.write(f"FAIL\t{x}\n")
    if fail:
        sys.exit("02_carrier_table.py: " + "; ".join(fail))
    with open(P + ".signature", "w") as f:
        f.write(sig + "\n")
    with open(P + ".DONE", "w") as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} people {len(uset_pat)} keys {nk}\n")
    print(f"02 {a.uset}: {len(uset_pat)} people, {nk} distinct keys, {sum(n_unique.values())} per-person records")


if __name__ == "__main__":
    main()
