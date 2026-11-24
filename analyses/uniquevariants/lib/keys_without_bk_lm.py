#!/usr/bin/env python3
"""A unique-set key table (extract_unique.py) without its BK_lm column, so that 04_classify_count.py reads it as a prepared key
table (it would take BK_lm for a region). The rows and every other column stay as they are. Used by xplat_jgt.sh.
Usage: python3 keys_without_bk_lm.py <in.keys.tsv.gz> <out.keys.tsv.gz>
"""
import os
import sys

import pandas as pd

src, dst = sys.argv[1], sys.argv[2]
df = pd.read_csv(src, sep="\t", dtype=str, keep_default_na=False)
if "BK_lm" not in df.columns:
    sys.exit(f"{src}: no BK_lm column (not a unique-set key table?)")
out = df.drop(columns=["BK_lm"])
os.makedirs(os.path.dirname(dst), exist_ok=True)
tmp = dst + ".tmp.gz"
out.to_csv(tmp, sep="\t", index=False, compression="gzip")
back = pd.read_csv(tmp, sep="\t", dtype=str, keep_default_na=False)
if len(back) != len(df) or list(back.columns) != list(out.columns):
    sys.exit(f"{dst}: written table differs ({len(back)} rows, {len(df)} expected)")
os.replace(tmp, dst)
print(f"{os.path.basename(dst)}: {len(back)} rows, BK_lm dropped")
