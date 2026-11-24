#!/usr/bin/env python3
"""The rule for region names (a plain identifier, unique, not the name of a record column), shared by lib/make_keys.py and
04_classify_count.py, and the prefix "stratum:" that 04 gives to the region columns of a key table.
"""
import re

COL_PREFIX = "stratum:"
PSEUDO = ("genome", "autosomes")
KEY_TABLE_COLS = ("CHROM", "POS", "REF", "ALT", "GT", "TYPE", "SPANDEL", "AUTO", "OLD_REC")
RECORD_COLS = ("key", "gkey", "side", "class", "cls", "BD", "BK", "BS", "line", "LMFLAG", "LM_FLAG", "BK_lm", "NA_MISSING",
               "SPANDEL_PAIR", "key_match", "vcfeval")
RESERVED = frozenset(PSEUDO + KEY_TABLE_COLS + RECORD_COLS)
_NAME = re.compile(r"^[A-Za-z0-9_]+$")


def check(names, where="strata"):
    """stop with a message when a stratum name is not a plain identifier, is repeated, or is a reserved column name"""
    names = list(names)
    bad = [n for n in names if not _NAME.match(str(n))]
    dup = sorted({n for n in names if names.count(n) > 1})
    clash = sorted({n for n in names if n in RESERVED})
    msg = []
    if bad:
        msg.append(f"names that are not plain identifiers: {bad}")
    if dup:
        msg.append(f"repeated names: {dup}")
    if clash:
        msg.append(f"names reserved for record columns or pseudo-strata: {clash}")
    if msg:
        raise SystemExit(f"{where}: " + "; ".join(msg))
    return names


def namespaced(df, names):
    """the key table with every stratum column renamed to COL_PREFIX + name"""
    return df.rename(columns={n: COL_PREFIX + n for n in names})
