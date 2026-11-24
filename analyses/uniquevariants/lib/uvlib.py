#!/usr/bin/env python3
"""uvlib.py - helpers shared by the Python scripts of this folder (imported, not run): the record key, the strata of a key table,
the counts per region and type (as in analyses/concordance/04_classify_count.py, which is loaded as a module) and the TSV readers.
"""
import csv
import functools
import gzip
import importlib.util
import os
import subprocess

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
UV_DIR = os.path.dirname(HERE)
CONC_DIR = os.environ.get("CONC_DIR", os.path.join(os.path.dirname(UV_DIR), "concordance"))
FLAG_COLS = ("AUTO", "SPANDEL", "BK_lm", "private_full", "private_full_excl_related", "ac1_unique", "ac1_unique_excl_related",
             "in_unique_excl_related")
TYPES = ["ALL", "SNV", "INDEL", "MNP", "OTHER"]
EX_MISSING = 64          # an input is missing, nothing written
EX_STALE = 65            # a DONE file with another fingerprint, nothing written
# strata whose rows in the concordance counts.tsv are not compared (none by default)
CONC_SKIP_STRATA = [s for s in os.environ.get("CONC_SKIP_STRATA", "").split(",") if s]


def conc04():
    """analyses/concordance/04_classify_count.py as a module (trim, side_records, classify)"""
    path = os.path.join(CONC_DIR, "04_classify_count.py")
    spec = importlib.util.spec_from_file_location("conc04", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_tsv(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        return list(csv.DictReader((ln for ln in f if not ln.startswith("#")), delimiter="\t"))


def tsv_row(path, keycol, key):
    for r in read_tsv(path):
        if r[keycol] == key:
            return r
    raise KeyError(f"{key} not in {path}")


def make_key(df):
    return df["CHROM"].astype(str) + ":" + df["POS"].astype(str) + ":" + df["REF"] + ":" + df["ALT"]


@functools.lru_cache(maxsize=None)
def strata_names():
    """the strata of the concordance config (STRATA), in its order; a stratum may not share its name with a flag column"""
    path = os.environ.get("STRATA", os.path.join(CONC_DIR, "config", "strata.tsv"))
    names = tuple(r["stratum"] for r in read_tsv(path))
    spec = importlib.util.spec_from_file_location("conc_strata_names", os.path.join(CONC_DIR, "lib", "strata_names.py"))
    rule = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rule)
    rule.check(names, f"uvlib: {path}")
    clash = sorted(set(names) & set(FLAG_COLS))
    if clash:
        raise SystemExit(f"uvlib: {path}: strata named like a flag column of this analysis: {clash}")
    return names


def load_keys(path, usecols=None):
    """a key table (prep/<set>_<platform>/<patient>.keys.tsv.gz or of this folder) with the 0/1 columns as booleans"""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, usecols=usecols)
    df["key"] = make_key(df)
    for c in df.columns:
        if c in FLAG_COLS or c in strata_names():
            df[c] = df[c] == "1"
    return df


def strata_columns(df):
    return [s for s in strata_names() if s in df.columns]


def stratum_mask(df, st):
    if st == "genome":
        return pd.Series(True, index=df.index)
    if st == "autosomes":
        return df["AUTO"].astype(bool)
    return df[st].astype(bool)


def count_by(df, strata, extra=None, value_cols=None):
    """rows stratum x type with n (and the sum of each boolean column of value_cols); extra = constant columns first"""
    rows = []
    for st in ["genome", "autosomes"] + list(strata):
        m = stratum_mask(df, st)
        for t in TYPES:
            mm = m if t == "ALL" else m & (df["TYPE"] == t)
            row = dict(extra or {})
            row.update({"stratum": st, "type": t, "n": int(mm.sum())})
            for c in value_cols or []:
                row[c] = int(df.loc[mm, c].astype(bool).sum())
            rows.append(row)
    return pd.DataFrame(rows)


def class_counts(df, strata, classcol, classes, extra=None):
    """rows stratum x type: n and one column per class (count of rows with that class); asserts the classes add up"""
    rows = []
    for st in ["genome", "autosomes"] + list(strata):
        m = stratum_mask(df, st)
        for t in TYPES:
            mm = m if t == "ALL" else m & (df["TYPE"] == t)
            sub = df.loc[mm, classcol].value_counts()
            row = dict(extra or {})
            row.update({"stratum": st, "type": t, "n": int(mm.sum())})
            for c in classes:
                row[c] = int(sub.get(c, 0))
            other = int(mm.sum()) - sum(row[c] for c in classes)
            if other:
                raise SystemExit(f"class_counts: {other} rows with a class outside {classes} ({st}, {t})")
            rows.append(row)
    return pd.DataFrame(rows)


def bcftools_lines(bcftools, args):
    """the stdout lines of a bcftools command (text), streamed"""
    p = subprocess.Popen([bcftools] + args, stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        yield line
    if p.wait() != 0:
        raise SystemExit(f"bcftools {' '.join(args)} failed")


def filter_vcf_by_keys(bcftools, vcf, keys, out_vcf_gz):
    """write the records of a prepared VCF whose CHROM:POS:REF:ALT is in keys (one ALT per record), bgzipped; returns the count"""
    n = 0
    with open(out_vcf_gz, "wb") as fo:
        bg = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=fo, text=True)
        for line in bcftools_lines(bcftools, ["view", vcf]):
            if line.startswith("#"):
                bg.stdin.write(line)
                continue
            f = line.split("\t", 5)
            if f"{f[0]}:{f[1]}:{f[3]}:{f[4]}" in keys:
                bg.stdin.write(line)
                n += 1
        bg.stdin.close()
        if bg.wait() != 0:
            raise SystemExit("bgzip failed")
    subprocess.run(["tabix", "-f", "-p", "vcf", out_vcf_gz], check=True)
    return n


def write_tsv(df, path, **kw):
    tmp = path + ".tmp"
    # mtime 0 in the gzip header: the same table gives the same bytes (md5) on a rerun
    comp = {"method": "gzip", "mtime": 0} if path.endswith(".gz") else None
    df.to_csv(tmp, sep="\t", index=False, na_rep="NA", compression=comp, **kw)
    os.replace(tmp, path)
