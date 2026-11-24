#!/usr/bin/env python3
"""The key table of one prepared call set: one row per record with its type (SNV, MNP, INDEL, OTHER), the SPANDEL flag (the record
held '*' before the split) and one 0/1 column per region (bcftools -T <BED> --targets-overlap variant). Prints counts on stdout.
Usage: python3 lib/make_keys.py --vcf prep.vcf.gz --strata strata.tsv --out prep.keys.tsv.gz [--bcftools bcftools]
"""
import argparse
import bisect
import io
import subprocess
import sys

import pandas as pd

import strata_names        # same folder (python3 lib/make_keys.py puts lib/ on the path)

AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
OK_BASES = set("ACGTNacgtn")


def run(cmd):
    """Run a command and return its stdout as text; stop on an error."""
    res = subprocess.run(cmd, check=True, stdout=subprocess.PIPE)
    return res.stdout.decode()


def var_type(ref, alt):
    """SNV, MNP, INDEL (any length change, complex alleles included) or OTHER (symbolic, '*', missing)."""
    if not ref or not alt or not set(ref) <= OK_BASES or not set(alt) <= OK_BASES:
        return "OTHER"
    if len(ref) == len(alt):
        return "SNV" if len(ref) == 1 else "MNP"
    return "INDEL"


def trim(pos, ref, alt):
    """minimal allele pair: common suffix, then common prefix, keeping one base on each side (as 04_classify_count.py)"""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return pos, ref, alt


def norm_shifts(df):
    """(left-shifted indels, those whose shift span holds the POS of another record of the same file)"""
    by_chrom = {}
    for c, p in zip(df["CHROM"], df["POS"]):
        by_chrom.setdefault(c, []).append(int(p))
    for v in by_chrom.values():
        v.sort()
    n_shift = n_cross = 0
    for c, p, o in zip(df["CHROM"], df["POS"], df["OLD_REC"]):
        f = o.split("|")
        if len(f) < 4 or not f[1].isdigit():
            continue
        alts = f[3].split(",")
        ia = int(f[4]) if len(f) > 4 and f[4].isdigit() else 1
        if not 1 <= ia <= len(alts):
            continue
        alt = alts[ia - 1]
        if not set(f[2]) <= OK_BASES or not alt or not set(alt) <= OK_BASES or len(f[2]) == len(alt):
            continue
        op, oref, _ = trim(int(f[1]), f[2], alt)
        p = int(p)
        if op <= p:
            continue
        n_shift += 1
        end = op + len(oref) - 1
        v = by_chrom[c]
        if bisect.bisect_right(v, end) - bisect.bisect_right(v, p) > 0:
            n_cross += 1
    return n_shift, n_cross


def gt_form(gt):
    al = gt.replace("|", "/").split("/")
    if "." in al:
        return "half_or_missing"
    if len(al) == 1:
        return "haploid"
    if len(al) == 2:
        if al[0] == al[1]:
            return "hom_alt" if al[0] != "0" else "hom_ref"
        return "het" if "0" in al else "het_alt_alt"
    return "polyploid"


def read_records(text, names):
    if not text:
        return pd.DataFrame({c: pd.Series(dtype=str) for c in names})
    return pd.read_csv(io.StringIO(text), sep="\t", header=None, names=names, dtype=str,
                       keep_default_na=False, quoting=3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--strata", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--bcftools", default="bcftools")
    a = ap.parse_args()

    df = read_records(run([a.bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%INFO/OLD_REC\n", a.vcf]),
                      ["CHROM", "POS", "REF", "ALT", "GT", "OLD_REC"])
    key = df["CHROM"] + ":" + df["POS"] + ":" + df["REF"] + ":" + df["ALT"]
    if key.duplicated().any():
        sys.exit(f"make_keys.py: {int(key.duplicated().sum())} duplicated keys in {a.vcf} (norm -d exact should remove them)")
    df["TYPE"] = [var_type(r, x) for r, x in zip(df["REF"], df["ALT"])]
    # SPANDEL: the record comes from a multi-allelic site that also held '*' (the allele next to an upstream deletion)
    # (OLD_REC is CHROM|POS|REF|ALT of the record before norm; only an allele can hold '*')
    df["SPANDEL"] = [1 if "*" in o else 0 for o in df["OLD_REC"]]
    df["AUTO"] = df["CHROM"].isin(AUTOSOMES).astype(int)

    strata = pd.read_csv(a.strata, sep="\t", header=None, names=["stratum", "bed"], dtype=str)
    strata_names.check(strata["stratum"], f"make_keys.py: {a.strata}")      # no stratum may share a name with a record column
    kidx = pd.Index(key)
    for name, bed in zip(strata["stratum"], strata["bed"]):
        s = read_records(run([a.bcftools, "query", "-T", bed, "--targets-overlap", "variant",
                              "-f", "%CHROM\t%POS\t%REF\t%ALT\n", a.vcf]), ["CHROM", "POS", "REF", "ALT"])
        skey = s["CHROM"] + ":" + s["POS"] + ":" + s["REF"] + ":" + s["ALT"]
        if not skey.isin(kidx).all():
            sys.exit(f"make_keys.py: stratum {name}: records not found in the main table")
        df[name] = key.isin(pd.Index(skey)).astype(int).values

    out_cols = ["CHROM", "POS", "REF", "ALT", "GT", "TYPE", "SPANDEL", "AUTO"] + list(strata["stratum"])
    df[out_cols].to_csv(a.out, sep="\t", index=False, compression="gzip")

    print(f"final_records\t{len(df)}")
    for t in ["SNV", "INDEL", "MNP", "OTHER"]:
        print(f"final_{t}\t{int((df['TYPE'] == t).sum())}")
    forms = pd.Series([gt_form(g) for g in df["GT"]], dtype=str)
    for f in ["het", "hom_alt", "half_or_missing", "haploid", "hom_ref", "het_alt_alt", "polyploid"]:
        print(f"gt_{f}\t{int((forms == f).sum())}")
    print(f"spandel_records\t{int(df['SPANDEL'].sum())}")
    n_shift, n_cross = norm_shifts(df)
    print(f"norm_left_shifted_indels\t{n_shift}")
    print(f"norm_shift_crosses_record\t{n_cross}")
    for name in strata["stratum"]:
        print(f"in_stratum_{name}\t{int(df[name].sum())}")


if __name__ == "__main__":
    main()
