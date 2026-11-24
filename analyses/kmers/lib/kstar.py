#!/usr/bin/env python3
"""kstar.py - k* analysis of the HG002 calls (kstar_hg002.sh): the copy number of the k-mers of every record in the reads and in
GRCh38, computed as merfin -dump does, the paralog-like flags (kstar, multicopy_ref, alt_elsewhere) and the Merfin x GIAB tables
split by them. kstar_hg002.sh uses the modes peak, windows, kmers, tables and verify (the copy numbers of the supported GIAB true
positives); compare and state were checks of the run.
"""
import argparse
import gzip
import os
import re
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "uniquevariants", "lib"))
import uvlib  # noqa: E402
from merfin_status import sample_states  # noqa: E402

GIAB_CODES = {"TP": 0, "FP": 1, "FP_gt": 2, "NE": 3}
GIAB_NAMES = {v: k for k, v in GIAB_CODES.items()}
FLAGS = ["kstar", "multicopy_ref", "alt_elsewhere", "paralog_like"]
BINS = [(-np.inf, -0.5, "lt_-0.5"), (-0.5, 0.5, "-0.5_to_0.5"), (0.5, 1.0, "0.5_to_1"), (1.0, np.inf, "ge_1")]
SENS_T = (0.5, 1.0)
SENS_METRICS = ("kstar_alt", "kstar_site", "either", "kstar_alt_raw", "kstar_site_raw", "either_raw")


def write(df, path, **kw):
    uvlib.write_tsv(df, path, float_format="%.6g", **kw)


# ---- peak ----
def read_hist(path):
    h = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2 and f[0].isdigit():
                h[int(f[0])] = int(f[1])
    return h


def hist_peak(h, fmax=5000):
    """the haploid peak from a k-mer histogram: after the first rise (the end of the error k-mers), the highest point of the
    3-bin running mean is the homozygous (diploid) peak, unless a peak at least half as high sits at 1.8-2.2 x of it"""
    if not h:
        return None
    xs = list(range(min(h), min(max(h), fmax) + 1))
    y = {x: h.get(x, 0) for x in xs}
    s = {x: (y.get(x - 1, y[x]) + y[x] + y.get(x + 1, y[x])) / 3.0 for x in xs}
    trough = None
    for x in xs[:-1]:
        if s[x + 1] > s[x]:
            trough = x
            break
    if trough is None:
        return None
    cand = [x for x in xs if x >= trough]
    m = max(cand, key=lambda x: s[x])
    near2 = [x for x in cand if 1.8 * m <= x <= 2.2 * m]
    if near2 and max(s[x] for x in near2) >= 0.5 * s[m]:
        d = max(near2, key=lambda x: s[x])
        return {"trough": trough, "highest": m, "diploid": float(d), "haploid": float(m),
                "rule": "the highest peak is the haploid one (a peak at 2 x of it)"}
    return {"trough": trough, "highest": m, "diploid": float(m), "haploid": m / 2.0,
            "rule": "the highest peak is the homozygous one; haploid = highest / 2"}


def gs_kcov(gsdir):
    p = os.path.join(gsdir, "model.txt")
    if not gsdir or not os.path.exists(p):
        return None
    with open(p) as fh:
        for line in fh:
            m = re.match(r"^kmercov\s+([0-9.eE+-]+)", line)
            if m:
                return float(m.group(1))
    return None


def lookup_ok(path):
    """the GenomeScope 2.0 lookup table (lines 'copies,probability'); recorded only, never given to merfin"""
    if not os.path.exists(path):
        return False
    n = 0
    with open(path) as fh:
        for line in fh:
            if not re.match(r"^\d+,[0-9.eE+-]+\s*$", line):
                return False
            n += 1
    return n > 0


def mode_peak(a):
    h = read_hist(a.hist)
    hp = hist_peak(h)
    if hp is None:
        sys.exit(f"kstar.py peak: no peak in {a.hist}")
    kc = gs_kcov(a.genomescope) if a.genomescope else None
    lk = os.path.join(a.genomescope, "lookup_table.txt") if a.genomescope else ""
    ratio = kc / hp["haploid"] if kc else float("nan")
    if kc and 0.75 <= ratio <= 1.33:
        method, peak = "genomescope", kc
    else:
        method, peak = "histogram", hp["haploid"]
    rows = [("haploid_peak", f"{peak:.2f}"), ("method", method), ("prob", "not used"),
            ("hist_trough", hp["trough"]), ("hist_highest", hp["highest"]), ("hist_diploid_peak", hp["diploid"]),
            ("hist_haploid_peak", f"{hp['haploid']:.2f}"), ("hist_rule", hp["rule"]),
            ("genomescope_kcov", f"{kc:.3f}" if kc else "NA"), ("genomescope_vs_hist", f"{ratio:.3f}" if kc else "NA"),
            ("genomescope_lookup_table", "valid" if a.genomescope and lookup_ok(lk) else "absent"),
            ("distinct_kmers", sum(h.values())), ("bins", len(h))]
    tmp = a.out + ".tmp"
    with open(tmp, "w") as fh:
        fh.writelines(f"{k}\t{v}\n" for k, v in rows)
    os.replace(tmp, a.out)
    print(f"{peak:.2f}\t{method}")


# ---- windows ----
def read_fai(path):
    out = {}
    with open(path + ".fai") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    return out


def read_contig(fa, ent):
    length, offset, lbases, lwidth = ent
    nlines = (length + lbases - 1) // lbases
    with open(fa, "rb") as fh:
        fh.seek(offset)
        raw = fh.read(nlines * lwidth)
    return raw.decode().replace("\n", "").replace("\r", "")[:length].upper()


def allele_windows(seq, pos, ref, alt, k):
    """(reference window, ALT window): the changed bases (REF and ALT without their shared leading, then trailing bases) with
    k - 1 reference bases on each side"""
    p0 = pos - 1
    if seq[p0:p0 + len(ref)] != ref.upper():
        raise SystemExit(f"kstar.py windows: REF of {pos}:{ref}:{alt} differs from the reference")
    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    r, x = ref[i:], alt[i:]
    j = 0
    while j < len(r) and j < len(x) and r[-1 - j] == x[-1 - j]:
        j += 1
    if j:
        r, x = r[:len(r) - j], x[:len(x) - j]
    s = p0 + i
    e = s + len(r)
    left, right = seq[max(0, s - (k - 1)):s], seq[e:e + (k - 1)]
    return left + r + right, left + x + right


def valid_kmers(w, k):
    if len(w) < k:
        return 0
    if not w.strip("ACGT"):
        return len(w) - k + 1
    bad = np.cumsum([0] + [c not in "ACGT" for c in w])
    return int(sum(1 for j in range(len(w) - k + 1) if bad[j + k] == bad[j]))


def mode_windows(a):
    fai = read_fai(a.ref)
    order = {c: i for i, c in enumerate(fai)}
    parts = []
    for p in a.keys.split(","):
        parts.append(pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, usecols=["CHROM", "POS", "REF", "ALT"]))
    df = pd.concat(parts, ignore_index=True).drop_duplicates()
    unknown = sorted(set(df["CHROM"]) - set(order))
    if unknown:
        sys.exit(f"kstar.py windows: contigs not in {a.ref}.fai: {unknown[:5]}")
    df["POSi"] = df["POS"].astype(int)
    df["o"] = df["CHROM"].map(order)
    df = df.sort_values(["o", "POSi", "REF", "ALT"]).reset_index(drop=True)
    idx_rows = []
    tmp = a.fasta + ".tmp"
    with open(tmp, "w") as fa:
        i = 0
        for chrom, g in df.groupby("o", sort=True):
            name = g["CHROM"].iloc[0]
            seq = read_contig(a.ref, fai[name])
            for pos, ref, alt in zip(g["POSi"], g["REF"], g["ALT"]):
                rw, aw = allele_windows(seq, pos, ref, alt, a.k)
                fa.write(f">{i}_R\n{rw}\n>{i}_A\n{aw}\n")
                idx_rows.append((i, f"{name}:{pos}:{ref}:{alt}", valid_kmers(rw, a.k), valid_kmers(aw, a.k)))
                i += 1
    os.replace(tmp, a.fasta)
    ix = pd.DataFrame(idx_rows, columns=["idx", "key", "ref_kmers", "alt_kmers"])
    write(ix, a.index)
    print(f"records\t{len(ix)}\tkmers_ref\t{int(ix['ref_kmers'].sum())}\tkmers_alt\t{int(ix['alt_kmers'].sum())}")


# ---- kmers: the counts of every window k-mer (meryl-lookup -wig-count) turned into the values of merfin -dump ----
WIG_HEAD = "variableStep chrom="
POS_W = 1 << 24      # row key = window ordinal x POS_W + position (a window is far shorter than 16 M bases)


def merfin_readk(count, peak):
    """readK of merfin-globals.C getK: 0 for a count of 0, 1 below the peak, else round(count / peak), half away from zero (C round:
    x - floor(x) is exact, so this equals it for every count)"""
    c = np.asarray(count, float)
    x = c / peak
    f = np.floor(x)
    return np.where(c == 0, 0.0, np.where(c < peak, 1.0, f + (x - f >= 0.5)))


def window_ordinal(names):
    """the ordinal of windows named <i>_R and <i>_A (windows.fa): 2 i for R, 2 i + 1 for A; -1 for any other name"""
    s = pd.Series(names, dtype=str)
    side, num = s.str[-2:], s.str[:-2]
    ok = (side.isin(["_R", "_A"]) & num.str.fullmatch(r"\d+")).to_numpy()
    o = np.full(len(s), -1, np.int64)
    o[ok] = num[ok].astype(np.int64).to_numpy() * 2 + (side[ok] == "_A").to_numpy()
    return o


def wig_chunks(path, chunk, info):
    """(key, count) arrays per chunk of a meryl-lookup -wig-count file of the windows: a header 'variableStep chrom=<window>' per
    window in the order of windows.fa, then '<1-based k-mer start>\t<count>' for each k-mer with a count above 0; key = window
    ordinal x POS_W + the 0-based start (merfin's seqPos). Stops unless every window comes once, in order, from 0, with increasing
    positions and whole positive counts; info[path] = the number of windows"""
    nxt, cur, last = 0, -1, -1
    for df in pd.read_csv(path, sep="\t", header=None, names=["a", "c"], dtype={"a": str, "c": float}, chunksize=chunk):
        hdr = df["c"].isna().to_numpy()
        a = df["a"].to_numpy()
        h = pd.Series(a[hdr], dtype=str)
        if not h.str.startswith(WIG_HEAD).all():
            sys.exit(f"kstar.py: {path}: a line without a count is not a '{WIG_HEAD}<window>' header")
        ords = window_ordinal(h.str.slice(len(WIG_HEAD)).to_numpy())
        if not np.array_equal(ords, np.arange(nxt, nxt + len(ords))):
            sys.exit(f"kstar.py: {path}: the windows are not those of windows.fa once each and in order (from ordinal {nxt})")
        nxt += len(ords)
        o = np.concatenate([[cur], ords])[np.cumsum(hdr)][~hdr]
        if len(ords):
            cur = int(ords[-1])
        if not len(o):
            continue
        pos = a[~hdr].astype(np.int64) - 1
        v = df["c"].to_numpy()[~hdr]
        if o[0] < 0 or (pos < 0).any() or (pos >= POS_W).any() or (v < 1).any() or (v != np.floor(v)).any():
            sys.exit(f"kstar.py: {path}: a count row before the first header, a position out of range or a count that is not a whole "
                     "number above 0")
        key = o * POS_W + pos
        if key[0] <= last or (np.diff(key) <= 0).any():
            sys.exit(f"kstar.py: {path}: positions not increasing within a window")
        last = int(key[-1])
        yield key, v.astype(np.int64)
    info[path] = nxt


def lookup_rows(reads_wig, ref_wig, chunk, info):
    """DataFrames (o, pos, count, asm) of every window k-mer with a count above 0 in the reads or in GRCh38 (the rows merfin -dump
    prints: readK, asmK or k* not 0), from the -wig-count files of the reduced read and GRCh38 databases, sorted by window ordinal o
    and position; each holds whole windows"""
    its = [wig_chunks(reads_wig, chunk, info), wig_chunks(ref_wig, chunk, info)]
    buf = [[np.empty(0, np.int64), np.empty(0, np.int64)] for _ in its]
    done = [False, False]
    while True:
        for i, it in enumerate(its):
            while not done[i] and (len(buf[i][0]) < chunk or buf[i][0][0] // POS_W == buf[i][0][-1] // POS_W):
                nx = next(it, None)
                if nx is None:
                    done[i] = True
                else:
                    buf[i] = [np.concatenate([buf[i][0], nx[0]]), np.concatenate([buf[i][1], nx[1]])]
        # windows below the last window of every buffer still being read are complete in both
        lim = [buf[i][0][-1] // POS_W * POS_W for i in (0, 1) if not done[i]]
        take = [b[0] < min(lim) if lim else np.ones(len(b[0]), bool) for b in buf]
        (rk, rv), (gk, gv) = [(b[0][t], b[1][t]) for b, t in zip(buf, take)]
        buf = [[b[0][~t], b[1][~t]] for b, t in zip(buf, take)]
        keys = np.union1d(rk, gk)
        if len(keys):
            count = np.zeros(len(keys), np.int64)
            asm = np.zeros(len(keys), np.int64)
            count[np.searchsorted(keys, rk)] = rv
            asm[np.searchsorted(keys, gk)] = gv
            yield pd.DataFrame({"o": keys // POS_W, "pos": keys % POS_W, "count": count, "asm": asm})
        if all(done) and not any(len(b[0]) for b in buf):
            return


def merfin_lines(o, pos, count, asm, peak):
    """the lines merfin -dump -peak <peak> prints for these k-mers: <i>_<R|A>, seqPos, readK, asmK, k* (%.2f), with readK = getK of the
    count and k* = getKmetric(readK, asmK) (inf when asmK is 0)"""
    rk = merfin_readk(count, peak)
    am = np.asarray(asm, float)
    ks = kmetric(rk, am)
    return [f"{i // 2}_{'A' if i % 2 else 'R'}\t{p}\t{r:.2f}\t{q:.2f}\t{z:.2f}"
            for i, p, r, q, z in zip(np.asarray(o).tolist(), np.asarray(pos).tolist(), rk.tolist(), am.tolist(), ks.tolist())]


def derived_windows(chunks, peak, keep=None):
    """(ordinal, [line, ...]) per window of lookup_rows chunks, with the lines merfin -dump would print; keep: the ordinals to give
    (sorted array; None = all)"""
    for ch in chunks:
        o = ch["o"].to_numpy()
        if keep is not None:
            m = np.isin(o, keep)
            ch, o = ch[m], o[m]
        if not len(o):
            continue
        lines = merfin_lines(o, ch["pos"].to_numpy(), ch["count"].to_numpy(), ch["asm"].to_numpy(), peak)
        st = np.flatnonzero(np.r_[True, o[1:] != o[:-1]])
        for s, e in zip(st.tolist(), np.r_[st[1:], len(o)].tolist()):
            yield int(o[s]), lines[s:e]


def dump_windows(path, stats):
    """(ordinal, [line, ...]) per window of a merfin -dump file. Complete lines only: reading stops at a last line without its newline
    or at the end of a cut gzip stream (stats: lines, cut)"""
    stats.update(lines=0, cut=0)
    cur, buf = None, []
    with gzip.open(path, "rt") as fh:
        try:
            for line in fh:
                if not line.endswith("\n"):
                    stats["cut"] = 1
                    break
                name = line.split("\t", 1)[0]
                i, _, s = name.rpartition("_")
                if s not in ("R", "A") or not i.isdigit():
                    sys.exit(f"kstar.py: {path}: the sequence name {name} is not a window of windows.fa")
                o = 2 * int(i) + (s == "A")
                if o != cur:
                    if cur is not None:
                        if o < cur:
                            sys.exit(f"kstar.py: {path}: the windows are not in the order of windows.fa ({name})")
                        yield cur, buf
                    cur, buf = o, []
                buf.append(line[:-1])
                stats["lines"] += 1
        except EOFError:
            stats["cut"] = 1
    if cur is not None:
        yield cur, buf


def compare_windows(derived, dumped, upto=None):
    """counts of a walk over two (ordinal, lines) streams sorted by ordinal (lines rebuilt from the counts, lines of merfin):
    windows in both, windows equal, rows of each, rows equal, windows in one only; upto: only windows below this ordinal.
    differences = rows of either side without an equal row on the other"""
    st = {"windows": 0, "windows_equal": 0, "rows_dump": 0, "rows_lookup": 0, "rows_equal": 0, "only_dump_windows": 0,
          "only_lookup_windows": 0, "examples": []}

    def nxt(it):
        x = next(it, None)
        return None if x is None or (upto is not None and x[0] >= upto) else x

    d, m = nxt(derived), nxt(dumped)
    while d is not None or m is not None:
        if m is None or (d is not None and d[0] < m[0]):
            st["only_lookup_windows"] += 1
            st["rows_lookup"] += len(d[1])
            st["examples"].append(f"only from the counts: {d[1][0]}")
            d = nxt(derived)
        elif d is None or m[0] < d[0]:
            st["only_dump_windows"] += 1
            st["rows_dump"] += len(m[1])
            st["examples"].append(f"only in the dump: {m[1][0]}")
            m = nxt(dumped)
        else:
            st["windows"] += 1
            st["rows_dump"] += len(m[1])
            st["rows_lookup"] += len(d[1])
            if d[1] == m[1]:
                st["windows_equal"] += 1
                st["rows_equal"] += len(m[1])
            else:
                st["rows_equal"] += len(set(d[1]) & set(m[1]))
                st["examples"].append(f"counts [{next(x for x in d[1] if x not in m[1]) if set(d[1]) - set(m[1]) else '-'}] "
                                      f"merfin [{next(x for x in m[1] if x not in d[1]) if set(m[1]) - set(d[1]) else '-'}]")
            d, m = nxt(derived), nxt(dumped)
    st["differences"] = st["rows_dump"] + st["rows_lookup"] - 2 * st["rows_equal"]
    st["examples"] = st["examples"][:5]
    return st


def compare_text(st):
    return "\t".join(f"{k}\t{st[k]}" for k in ("windows", "windows_equal", "rows_dump", "rows_lookup", "rows_equal",
                                                 "only_dump_windows", "only_lookup_windows", "differences"))


def summarise(ch, peak):
    """per window ordinal of a lookup_rows chunk: rows, present (in the reads), one (single-copy), x, x_raw, readK (medians), a1, a2
    (k-mers 1 or more and 2 or more times in GRCh38)"""
    o = ch["o"].to_numpy()
    c = ch["count"].to_numpy(float)
    ak = ch["asm"].to_numpy(float)
    rk = merfin_readk(c, peak)
    isalt = (o % 2) == 1
    pres = rk > 0
    # single-copy k-mers: ALT absent from GRCh38 and seen in the reads; REF once in GRCh38 (an absent one counts as 0)
    one = np.where(isalt, (ak == 0) & pres, ak == 1)
    t = pd.DataFrame({"o": o, "present": pres, "one": one, "x": np.where(one, rk, np.nan), "x_raw": np.where(one, c / peak, np.nan),
                      "readK": np.where(pres, rk, np.nan), "a1": ak >= 1, "a2": ak >= 2, "rows": 1})
    return t.groupby("o", sort=False).agg(rows=("rows", "sum"), present=("present", "sum"), one=("one", "sum"), x=("x", "median"),
                                          x_raw=("x_raw", "median"), readK=("readK", "median"), a1=("a1", "sum"), a2=("a2", "sum"))


def fasta_windows(path):
    with open(path) as fh:
        o = window_ordinal([x[1:].split()[0] for x in fh if x.startswith(">")])
    if (o < 0).any():
        sys.exit(f"kstar.py: {path}: a sequence name is not a window of windows.fa")
    return np.unique(o)


def mode_kmers(a):
    ix = pd.read_csv(a.index, sep="\t", dtype={"idx": np.int64, "key": str, "ref_kmers": np.int64, "alt_kmers": np.int64})
    if not np.array_equal(ix["idx"].to_numpy(), np.arange(len(ix))):
        sys.exit(f"kstar.py kmers: {a.index} is not numbered 0, 1, 2, ...")
    peak = float(a.peak)
    spot = fasta_windows(a.spot_fasta)
    info, res, srows, nrows, t0 = {}, [], [], 0, time.time()
    for ch in lookup_rows(a.reads_wig, a.ref_wig, a.chunk, info):
        nrows += len(ch)
        res.append(summarise(ch, peak))
        srows.append(ch[np.isin(ch["o"].to_numpy(), spot)])
        print(f"progress\tk-mer rows {nrows}\twindows {int(ch['o'].iat[-1]) + 1} of {2 * len(ix)}\t{time.time() - t0:.0f} s", flush=True)
    for p in (a.reads_wig, a.ref_wig):
        if info.get(p) != 2 * len(ix):
            sys.exit(f"kstar.py kmers: {p} lists {info.get(p)} windows, {a.index} {2 * len(ix)}")
    g = pd.concat(res) if res else pd.DataFrame(columns=["rows", "present", "one", "x", "x_raw", "readK", "a1", "a2"])
    if g.index.duplicated().any():
        sys.exit("kstar.py kmers: the rows of one window are not together")
    g = g.reset_index()
    g["idx"], g["side"] = g["o"] // 2, np.where(g["o"] % 2 == 1, "A", "R")
    out = ix.copy()
    for s, pre in (("R", "ref"), ("A", "alt")):
        sub = g[g["side"] == s].set_index("idx")
        out[f"{pre}_rows"] = out["idx"].map(sub["rows"]).fillna(0).astype(int)
        out[f"{pre}_present"] = out["idx"].map(sub["present"]).fillna(0).astype(int)
        out["ref_once" if s == "R" else "alt_new"] = out["idx"].map(sub["one"]).fillna(0).astype(int)
        out[f"{pre}_readK"] = out["idx"].map(sub["readK"])
        out[f"{pre}_x"] = out["idx"].map(sub["x"])
        out[f"{pre}_x_raw"] = out["idx"].map(sub["x_raw"])
        out[f"{pre}_a1"] = out["idx"].map(sub["a1"]).fillna(0).astype(int)
        out[f"{pre}_{'multi' if s == 'R' else 'in_ref'}"] = out["idx"].map(sub["a2" if s == "R" else "a1"]).fillna(0).astype(int)
    over = (out["ref_rows"] > out["ref_kmers"]) | (out["alt_rows"] > out["alt_kmers"])
    if over.any():
        sys.exit(f"kstar.py kmers: {int(over.sum())} windows with more k-mer rows than valid k-mers (wrong index?)")
    short = out["ref_a1"] != out["ref_kmers"]
    if short.any():
        sys.exit(f"kstar.py kmers: {int(short.sum())} REF windows with k-mers not found in GRCh38 (every REF k-mer is in GRCh38; "
                 f"wrong index or GRCh38 database?), e.g. {out.loc[short, 'key'].head(3).tolist()}")
    sst = {}
    sp = pd.concat(srows, ignore_index=True) if srows else pd.DataFrame(columns=["o", "pos", "count", "asm"])
    cs = compare_windows(derived_windows(iter([sp]), peak), dump_windows(a.spot_dump, sst))
    if sst["cut"] or cs["differences"] or cs["only_dump_windows"]:
        sys.exit(f"kstar.py kmers: the spot check failed: the lines rebuilt from the counts differ from merfin -dump -peak {peak} "
                 f"({a.spot_dump}; cut {sst['cut']}): {compare_text(cs)}; {cs['examples']}")
    cols = ["key", "ref_kmers", "ref_present", "ref_once", "ref_readK", "ref_x", "ref_x_raw", "ref_multi", "alt_kmers",
            "alt_present", "alt_new", "alt_readK", "alt_x", "alt_x_raw", "alt_in_ref"]
    write(out[cols], a.out)
    print(f"spot check: {len(spot)} windows of {a.spot_fasta}, merfin -dump -peak {peak}: {compare_text(cs)}")
    print(f"records\t{len(out)}\tkmer_rows\t{nrows}\tno_rows_alt\t{int((out['alt_rows'] == 0).sum())}"
          f"\tref_kmers_in_grch38\t{len(out)}\tspot_windows\t{len(spot)}\tspot_rows\t{cs['rows_dump']}\tspot_differences\t{cs['differences']}")


def mode_compare(a):
    peak = float(a.peak)
    upto = None
    if a.partial:
        st0, last = {}, None
        for last, _ in dump_windows(a.dump, st0):
            pass
        if last is None:
            sys.exit(f"kstar.py compare: no complete line in {a.dump}")
        upto = last
        print(f"partial dump: {st0['lines']} complete lines (cut line or stream: {st0['cut']}); the last window, ordinal {last}, may be "
              f"cut: windows below it are compared")
    info, sd = {}, {}
    st = compare_windows(derived_windows(lookup_rows(a.reads_wig, a.ref_wig, a.chunk, info), peak), dump_windows(a.dump, sd), upto)
    print(f"lookup windows: reads {info.get(a.reads_wig)}, GRCh38 {info.get(a.ref_wig)}; dump lines read {sd['lines']}")
    print(f"compare\t{compare_text(st)}")
    for x in st["examples"]:
        print(f"  {x}")
    ok = st["differences"] == 0 and st["windows"] > 0 and info.get(a.reads_wig) == info.get(a.ref_wig)
    print("COMPARE_EQUAL" if ok else "COMPARE_DIFFERENT")
    sys.exit(0 if ok else 1)


# ---- tables ----
def dosage(gt):
    d = sum(1 for x in re.split(r"[/|]", gt) if x == "1")
    return d if d > 0 else 1


def kmetric(o, e):
    """k* of merfin (merfin-globals.H getKmetric): 0 when o = 0; o / e - 1 when o > e (inf when e = 0); -(e / o - 1) when o < e"""
    o = np.asarray(o, float)
    e = np.asarray(e, float)
    out = np.zeros(len(o))
    m = o > 0
    hi = m & (o > e)
    lo = m & (o < e)
    with np.errstate(divide="ignore"):
        out[hi] = o[hi] / e[hi] - 1
    out[lo] = -(e[lo] / o[lo] - 1)
    return out


def status(path):
    return pd.read_csv(path, sep="\t", dtype={"key": str, "merfin": int}, keep_default_na=False)


def add_kstar(df, haploid_xy, min_kmers=5):
    """site_copies c (2; 1 on chrY and chrX outside XPAR with haploid_xy), exp_alt = min(d, c), and kstar_alt, kstar_site and
    their _raw versions (from the unrounded medians) on a table with CHROM, d, ref_once, ref_x, ref_x_raw, alt_new, alt_x,
    alt_x_raw [, XPAR]; kstar_alt is NA when alt_new < min_kmers, kstar_site also when ref_once < min_kmers"""
    site = np.full(len(df), 2.0)
    if haploid_xy:
        isx = (df["CHROM"] == "chrX").to_numpy()
        if "XPAR" not in df.columns and isx.any():
            sys.exit("kstar.py tables: --haploid-xy needs the XPAR column of the key table for the chrX records")
        xpar = df["XPAR"].astype(bool).to_numpy() if "XPAR" in df.columns else np.zeros(len(df), bool)
        site[((df["CHROM"] == "chrY").to_numpy() | isx) & ~xpar] = 1.0
    df["site_copies"] = site
    df["exp_alt"] = np.minimum(df["d"].to_numpy(float), site)
    oka = (df["alt_new"] >= min_kmers).to_numpy()
    oks = oka & (df["ref_once"] >= min_kmers).to_numpy()
    for sfx in ("", "_raw"):
        xa = df[f"alt_x{sfx}"].fillna(0).clip(lower=0).to_numpy(float)
        xr = df[f"ref_x{sfx}"].fillna(0).clip(lower=0).to_numpy(float)
        df[f"kstar_alt{sfx}"] = np.where(oka, kmetric(xa, df["exp_alt"].to_numpy()), np.nan)
        df[f"kstar_site{sfx}"] = np.where(oks, kmetric(xr + xa, site), np.nan)
    return df


def full_records(keys_path, st_path, truth_path, ks, kmin, haploid_xy, min_kmers):
    k = uvlib.load_keys(keys_path, usecols=lambda c: c != "SPANDEL").drop(columns=["POS", "REF", "ALT"])
    df = k.merge(status(st_path), on="key", how="left")
    if df["merfin"].isna().any():
        sys.exit(f"kstar.py tables: {int(df['merfin'].isna().sum())} keys of {keys_path} have no merfin status")
    tr = pd.read_csv(truth_path, sep="\t", dtype=str, keep_default_na=False, usecols=["key", "status"])
    df = df.merge(tr, on="key", how="left")
    df["giab"] = df["status"].fillna("NE")
    bad = sorted(set(df["giab"]) - set(GIAB_CODES))
    if bad:
        sys.exit(f"kstar.py tables: unknown GIAB status {bad} in {truth_path}")
    df = df.merge(ks, on="key", how="left")
    if df["alt_kmers"].isna().any():
        sys.exit(f"kstar.py tables: {int(df['alt_kmers'].isna().sum())} keys of {keys_path} are not in the k-mer table")
    df["d"] = [dosage(g) for g in df["GT"]]
    add_kstar(df, haploid_xy, min_kmers)
    add_flags(df, kmin)
    df["merfin"] = df["merfin"].astype(int)
    return df


def add_flags(df, kmin):
    """the flags kstar (NA with kstar_alt), multicopy_ref, alt_elsewhere, paralog_like on a table with kstar_alt, kstar_site,
    ref_kmers, ref_multi, alt_kmers, alt_in_ref, alt_present, alt_new. alt_elsewhere: every ALT k-mer is in GRCh38, or every ALT
    k-mer seen in the reads is (a cluster of paralog differences: the k-mers of one call alone are those of the other copy, the
    k-mers that also hold the REF base of a neighbour are neither in GRCh38 nor in the reads)"""
    hit = (df["kstar_alt"] >= kmin) | (df["kstar_site"] >= kmin)
    df["kstar"] = np.where(hit, 1.0, np.where(df["kstar_alt"].isna(), np.nan, 0.0))
    df["multicopy_ref"] = ((df["ref_kmers"] > 0) & (df["ref_multi"] > df["ref_kmers"] / 2)).astype(int)
    df["alt_elsewhere"] = (((df["alt_kmers"] > 0) & (df["alt_in_ref"] == df["alt_kmers"]))
                           | ((df["alt_present"] > 0) & (df["alt_new"] == 0))).astype(int)
    df["paralog_like"] = (hit | (df["multicopy_ref"] == 1) | (df["alt_elsewhere"] == 1)).astype(int)
    return df


def slices(df):
    strata = ["genome", "autosomes"] + uvlib.strata_columns(df)
    for st in strata:
        m = uvlib.stratum_mask(df, st).to_numpy()
        for t in uvlib.TYPES:
            yield st, t, (m if t == "ALL" else m & (df["TYPE"] == t).to_numpy())


def sens_values(df):
    v = {m: df[m].to_numpy() for m in ("kstar_alt", "kstar_site", "kstar_alt_raw", "kstar_site_raw")}
    v["either"] = np.fmax(v["kstar_alt"], v["kstar_site"])
    v["either_raw"] = np.fmax(v["kstar_alt_raw"], v["kstar_site_raw"])
    return v


FLAG_VALUES = {0: "0", 1: "1", 2: "NA"}   # flag code 2 = NA (kstar only)
ANY = (0, 1, 2)


def ncount(c, fl, me, gs):
    """records of a bincount (code = flag x 8 + merfin x 4 + giab) with flag in fl, merfin in me, giab in gs"""
    return int(sum(c[f * 8 + m * 4 + GIAB_CODES[g]] for f in fl for m in me for g in gs))


def aggregate(df, ex, x2, summ, bins, sens):
    g = df["giab"].map(GIAB_CODES).to_numpy()
    mf = df["merfin"].to_numpy()
    ka = df["kstar_alt"].to_numpy()
    sv = sens_values(df)
    sfp = (mf == 1) & (g >= 1) & (g <= 2)
    stp = (mf == 1) & (g == 0)
    codes = {f: df[f].fillna(2).to_numpy().astype(int) * 8 + mf * 4 + g for f in FLAGS}
    ev_g, fp_g = ("TP", "FP", "FP_gt"), ("FP", "FP_gt")
    for st, t, m in slices(df):
        cnt = {f: np.bincount(codes[f][m], minlength=24) for f in FLAGS}
        for f in FLAGS:
            for c in np.nonzero(cnt[f])[0]:
                x2.append({**ex, "stratum": st, "type": t, "flag_name": f, "flag": FLAG_VALUES[int(c // 8)],
                           "merfin": int((c // 4) % 2), "giab": GIAB_NAMES[int(c % 4)], "n": int(cnt[f][c])})
        c = cnt["kstar"]
        rec = int(c.sum())
        ev = ncount(c, ANY, (0, 1), ev_g)
        tp = ncount(c, ANY, (0, 1), ("TP",))
        fp = ev - tp
        mev = ncount(c, ANY, (1,), ev_g)
        mfp = ncount(c, ANY, (1,), fp_g)
        mtp = ncount(c, ANY, (1,), ("TP",))
        prec = tp / ev if ev else np.nan
        row = {**ex, "stratum": st, "type": t, "records": rec, "evaluable": ev, "giab_tp": tp, "giab_precision": prec,
               "merfin_share_evaluable": mev / ev if ev else np.nan, "gap_merfin_minus_precision": (mev / ev - prec) if ev else np.nan,
               "giab_fp": fp, "giab_fp_merfin": mfp, "merfin_supports_giab_fp": mfp / fp if fp else np.nan}
        for f in ("kstar", "paralog_like"):
            cf = cnt[f]
            ffp = ncount(cf, (1,), (1,), fp_g)
            ftp = ncount(cf, (1,), (1,), ("TP",))
            unf = ncount(cf, (0, 2), (1,), ev_g)
            row.update({f"flagged_{f}": ncount(cf, (1,), (0, 1), tuple(GIAB_CODES)),
                        f"giab_fp_merfin_flagged_{f}": ffp,
                        f"share_giab_fp_merfin_flagged_{f}": ffp / mfp if mfp else np.nan,
                        f"giab_tp_merfin_flagged_{f}": ftp,
                        f"share_giab_tp_merfin_flagged_{f}": ftp / mtp if mtp else np.nan,
                        f"merfin_unflagged_share_evaluable_{f}": unf / ev if ev else np.nan,
                        f"gap_unflagged_{f}": (unf / ev - prec) if ev else np.nan})
        row.update({"kstar_na": ncount(c, (2,), (0, 1), tuple(GIAB_CODES)), "giab_fp_merfin_kstar_na": ncount(c, (2,), (1,), fp_g),
                    "giab_tp_merfin_kstar_na": ncount(c, (2,), (1,), ("TP",))})
        summ.append(row)
        nfp, ntp = int((m & sfp).sum()), int((m & stp).sum())
        for mt in SENS_METRICS:
            for th in SENS_T:
                fl = m & (sv[mt] >= th)
                a_fp, a_tp = int((fl & sfp).sum()), int((fl & stp).sum())
                sens.append({**ex, "stratum": st, "type": t, "metric": mt, "threshold": th, "flagged": int(fl.sum()),
                             "giab_fp_merfin": nfp, "giab_fp_merfin_flagged": a_fp, "share_fp_flagged": a_fp / nfp if nfp else np.nan,
                             "giab_tp_merfin": ntp, "giab_tp_merfin_flagged": a_tp, "share_tp_flagged": a_tp / ntp if ntp else np.nan})
        sup = m & (mf == 1)
        for lo, hi, lab in BINS + [(None, None, "NA")]:
            vc = np.bincount(g[sup & (np.isnan(ka) if lab == "NA" else (ka >= lo) & (ka < hi))], minlength=4)
            for gi in range(4):
                if vc[gi]:
                    bins.append({**ex, "stratum": st, "type": t, "kstar_alt_bin": lab, "giab": GIAB_NAMES[gi], "n": int(vc[gi])})


README = """Outputs of the k* analysis of HG002 (kstar_hg002.sh, lib/kstar.py; one folder per platform)

records_<set>.tsv.gz  one row per record of a full set: key, d (ALT dosage of GT), merfin (1 = kept by merfin -filter), giab (TP, FP,
                      FP_gt; NE = outside the v5.0q benchmark regions or not evaluable), ref_kmers / alt_kmers (valid k-mers of the
                      two windows), ref_present / alt_present (seen in the reads), ref_readK / alt_readK (median read copies of
                      the k-mers seen), ref_once (REF k-mers once in GRCh38), ref_x (their median readK, absent = 0: REF copies
                      at the site), alt_new (ALT k-mers absent from GRCh38 and seen in the reads), alt_x (their median readK: ALT
                      copies), ref_multi (REF k-mers with 2 or more GRCh38 occurrences), alt_in_ref (ALT k-mers found in GRCh38),
                      site_copies (expected copies of the site), exp_alt (expected ALT copies, min(d, site_copies)), kstar_alt
                      (NA when alt_new < 5), kstar_site (NA also when ref_once < 5), the unrounded ref_x_raw, alt_x_raw,
                      kstar_alt_raw, kstar_site_raw (extra, not used by the flags), and the flags kstar (NA with kstar_alt),
                      multicopy_ref, alt_elsewhere (all ALT k-mers in GRCh38, or all of those seen in the reads), paralog_like
                      (analyses/kmers/README.md, section "k* on HG002")
kstar_2x2.tsv         uset (all_<set> or <unique set>_<platform>), platform, mode (context, direct), stratum, type, flag_name, flag
                      (0, 1, NA), merfin, giab, n: the Merfin x GIAB 2x2 split by each flag (rows with n > 0)
kstar_summary.tsv     per uset, platform, mode, stratum, type: records; evaluable (TP + FP + FP_gt); giab_tp; giab_precision;
                      merfin_share_evaluable; gap_merfin_minus_precision; giab_fp; giab_fp_merfin (false positives Merfin keeps);
                      merfin_supports_giab_fp; and for the flags kstar and paralog_like: flagged_<f> (all records),
                      giab_fp_merfin_flagged_<f> and its share of giab_fp_merfin, giab_tp_merfin_flagged_<f> and its share of the
                      Merfin-supported TP (the cost of the flag), merfin_unflagged_share_evaluable_<f> and gap_unflagged_<f> (the
                      Merfin share and the gap when flagged support is not counted as support); kstar_na,
                      giab_fp_merfin_kstar_na, giab_tp_merfin_kstar_na (records without k*)
kstar_bins.tsv        the kstar_alt of the Merfin-supported records in 4 bins (lt_-0.5, -0.5_to_0.5, 0.5_to_1, ge_1) and NA, by GIAB
                      status
kstar_sensitivity.tsv uset, platform, mode, stratum, type, metric (kstar_alt, kstar_site, either = the kstar rule, and the _raw
                      versions), threshold (0.5, 1): flagged records, Merfin-supported GIAB FP and TP, how many of them are flagged
                      and the shares
peak.tsv              the haploid peak (readK = round(count / peak) as with merfin -peak) and its method (genomescope or
                      histogram); -prob is not used
KSTAR_DONE, verify.txt  the end of the run (kstar_hg002.sh) and the checks of kstar.py verify
"""


def mode_tables(a):
    sets = pd.read_csv(a.sets, sep="\t", dtype=str, keep_default_na=False).set_index("set_id")
    samples = pd.read_csv(os.path.join(a.sheet_dir, "samples.tsv"), sep="\t", dtype=str, keep_default_na=False)
    vcfs = pd.read_csv(os.path.join(a.sheet_dir, "vcfs.tsv"), sep="\t", dtype=str, keep_default_na=False)
    me = samples[(samples["platform"] == a.platform) & (samples["sample"] == a.sample)]
    if len(me) != 1:
        sys.exit(f"kstar.py tables: {a.platform} {a.sample} is not once in the sheet")
    state = sample_states(me, a.results, a.state)[(a.platform, a.sample)]
    if state != "current":
        print(f"MISSING_INPUT\tsample\t{a.platform}:{a.sample}\tstate {state} (run merfin_sample.sh on it first)")
        sys.exit(64)
    res = os.path.join(a.results, a.platform, a.sample)
    rows = vcfs[(vcfs["platform"] == a.platform) & (vcfs["sample"] == a.sample)]
    ks = pd.read_csv(a.kstats, sep="\t", dtype={"key": str}, keep_default_na=False, na_values=["NA"])
    os.makedirs(a.out, exist_ok=True)
    full = {}
    x2, summ, bins, sens = [], [], [], []
    for r in rows[rows["kind"] == "all"].itertuples(index=False):
        hname = sets.loc[r.set_id, "hg002_name"]
        tp = os.path.join(a.uv_out, "hg002", "truth", f"{hname}_{a.platform}.status.tsv.gz")
        if not os.path.exists(tp):
            print(f"MISSING_INPUT\ttruth\t{a.platform}:{a.sample}\t{tp}")
            sys.exit(64)
        df = full_records(r.keys, os.path.join(res, f"{r.set_id}.status.tsv.gz"), tp, ks, a.kstar_min, a.haploid_xy, a.min_kmers)
        full[r.set_id] = df
        keep = ["key", "d", "merfin", "giab", "ref_kmers", "ref_present", "ref_readK", "ref_once", "ref_x", "ref_multi",
                "alt_kmers", "alt_present", "alt_readK", "alt_new", "alt_x", "alt_in_ref", "site_copies", "exp_alt", "kstar_alt",
                "kstar_site", "ref_x_raw", "alt_x_raw", "kstar_alt_raw", "kstar_site_raw"] + FLAGS
        write(df[keep], os.path.join(a.out, f"records_{hname}.tsv.gz"))
        aggregate(df, {"uset": f"all_{hname}", "platform": a.platform, "mode": "context"}, x2, summ, bins, sens)
    for r in rows[rows["kind"] == "unique"].itertuples(index=False):
        base = full[r.context]
        keys = set(uvlib.load_keys(r.keys, usecols=["CHROM", "POS", "REF", "ALT"])["key"])
        sub = base[base["key"].isin(keys)].copy()
        if len(sub) != len(keys):
            sys.exit(f"kstar.py tables: {len(keys) - len(sub)} keys of {r.keys} are not in the full set {r.context}")
        uset = f"{sets.loc[r.set_id, 'hg002_name']}_{a.platform}"
        aggregate(sub, {"uset": uset, "platform": a.platform, "mode": "context"}, x2, summ, bins, sens)
        dp = os.path.join(res, f"{r.set_id}.status.tsv.gz")
        if r.run == "yes" and os.path.exists(dp):
            d = status(dp)
            if set(d["key"]) != keys:
                sys.exit(f"kstar.py tables: the direct status of {r.set_id} and its key table hold different keys")
            sub = sub.drop(columns=["merfin"]).merge(d, on="key", how="left")
            sub["merfin"] = sub["merfin"].astype(int)
            aggregate(sub, {"uset": uset, "platform": a.platform, "mode": "direct"}, x2, summ, bins, sens)
    write(pd.DataFrame(x2), os.path.join(a.out, "kstar_2x2.tsv"))
    write(pd.DataFrame(summ), os.path.join(a.out, "kstar_summary.tsv"))
    write(pd.DataFrame(bins), os.path.join(a.out, "kstar_bins.tsv"))
    write(pd.DataFrame(sens), os.path.join(a.out, "kstar_sensitivity.tsv"))
    with open(os.path.join(a.out, "README.txt"), "w") as fh:
        fh.write(README)
    for sid, df in full.items():
        sup = df[df["merfin"] == 1]
        fpm = sup[sup["giab"].isin(["FP", "FP_gt"])]
        tpm = sup[sup["giab"] == "TP"]
        share = (tpm["kstar"] == 1).mean() if len(tpm) else float("nan")
        print(f"{sid}\trecords\t{len(df)}\tmerfin\t{len(sup)}\tgiab_fp_merfin\t{len(fpm)}\tof_them_kstar\t{int((fpm['kstar'] == 1).sum())}"
              f"\tof_them_kstar_na\t{int(fpm['kstar'].isna().sum())}\tof_them_paralog_like\t{int(fpm['paralog_like'].sum())}"
              f"\tgiab_tp_merfin\t{len(tpm)}\tshare_tp_kstar\t{share:.4f}\tkstar_na\t{int(df['kstar'].isna().sum())}")


AUTOSOME = r"^(chr)?([1-9]|1[0-9]|2[0-2])$"
CALIB_X = {1: (0.85, 1.15), 2: (1.7, 2.3)}   # median alt_x_raw of true supported calls in single-copy context, by dosage


def isolated(keys, k=31):
    """True for a record with no other record of the table within k - 1 bases (its ALT window holds no other call; in a cluster
    every ALT k-mer of one call can cover a neighbour, so none is in the reads although merfin supports the combination)"""
    p = keys.str.split(":", expand=True)
    t = pd.DataFrame({"c": p[0], "s": p[1].astype(np.int64)})
    t["e"] = t["s"] + p[2].str.len() - 1
    t = t.sort_values(["c", "s"], kind="stable")
    prev_end = t.groupby("c", sort=False)["e"].transform(lambda x: x.cummax().shift())
    next_start = t.groupby("c", sort=False)["s"].shift(-1)
    iso = (prev_end.isna() | (prev_end < t["s"] - (k - 1))) & (next_start.isna() | (next_start > t["e"] + (k - 1)))
    return iso.reindex(keys.index).to_numpy()


def calib_mask(rec):
    """the calibration records of one records_<set>.tsv.gz table: Merfin-supported GIAB TP with multicopy_ref 0 and
    alt_elsewhere 0 on the autosomes, with no other record within 30 bases"""
    chrom = rec["key"].str.split(":").str[0]
    return ((rec["merfin"] == 1) & (rec["giab"] == "TP") & (rec["multicopy_ref"] == 0) & (rec["alt_elsewhere"] == 0)
            & chrom.str.match(AUTOSOME) & isolated(rec["key"])).to_numpy()


def calibration(rec, calib_min):
    """(ok, detail) of the copy calibration on the calib_mask records; a wrong peak moves the medians by its factor. The detail
    ends each failed part with [fails: ...]"""
    s = rec[calib_mask(rec)]
    ok, parts = True, []
    for d, (lo, hi) in CALIB_X.items():
        sd = s[s["d"] == d]
        mx, mk = sd["alt_x_raw"].median(), sd["alt_readK"].median()
        bad = [w for w, b in (("n", len(sd) < calib_min), ("alt_x_raw", not lo <= mx <= hi),
                              ("alt_readK", not abs(mk - d) <= 0.25)) if b]
        ok &= not bad
        parts.append(f"d {d}: n {len(sd)} (>= {calib_min}), median alt_x_raw {mx:.3f} ({lo}-{hi}), median alt_readK {mk:.2f}"
                     + (f" [fails: {', '.join(bad)}]" if bad else ""))
    absent = float((s["alt_present"] == 0).mean()) if len(s) else float("nan")
    ok &= bool(absent < 0.01)
    parts.append(f"alt_present 0: {absent:.4f} (< 0.01)" + ("" if absent < 0.01 else " [fails: alt_present]"))
    return bool(ok), "; ".join(parts)


def mode_verify(a):
    nfail = nwarn = 0

    def line(name, ok, detail, warn=False):
        nonlocal nfail, nwarn
        res = "PASS" if ok else ("WARN" if warn else "FAIL")
        print(f"VERIFY\t{name}\t{res}\t{detail}")
        nfail += res == "FAIL"
        nwarn += res == "WARN"

    x2 = pd.read_csv(os.path.join(a.kstar, "kstar_2x2.tsv"), sep="\t", dtype=str)
    x2["n"] = x2["n"].astype(int)
    k = x2[(x2["platform"] == a.platform) & (x2["flag_name"] == "kstar") & (x2["stratum"] == "genome") & (x2["type"] == "ALL")]
    k = k.groupby(["uset", "mode", "merfin", "giab"])["n"].sum().reset_index()
    m3 = pd.read_csv(os.path.join(a.hg002, "merfin_giab_2x2.tsv"), sep="\t", dtype=str)
    m3 = m3[m3["platform"] == a.platform].copy()
    m3["n"] = m3["n"].astype(int)
    # suffix _03 (the name of the run step): the counts of hg002/merfin_giab_2x2.tsv; _k: the counts of kstar_2x2.tsv
    j = m3.merge(k, on=["uset", "mode", "merfin", "giab"], how="outer", suffixes=("_03", "_k"))
    diff = j[j["n_03"].fillna(-1) != j["n_k"].fillna(-1)]
    line("2x2_equals_merfin_giab_2x2", len(m3) > 0 and len(diff) == 0,
         f"{a.platform}: {len(m3)} rows of hg002/merfin_giab_2x2.tsv, {len(k)} of kstar_2x2 (flag summed), {len(diff)} differ"
         + (f": {diff.head(3).to_dict('records')}" if len(diff) else ""))
    sm = pd.read_csv(os.path.join(a.kstar, "kstar_summary.tsv"), sep="\t")
    se = pd.read_csv(os.path.join(a.kstar, "kstar_sensitivity.tsv"), sep="\t")
    sel = lambda d: d[(d["platform"] == a.platform) & d["uset"].str.startswith("all_") & (d["mode"] == "context")  # noqa: E731
                      & (d["stratum"] == "genome") & (d["type"] == "ALL")]
    g, gs = sel(sm), sel(se)
    line("full_sets", len(g) == 3, f"{len(g)} full sets in kstar_summary.tsv (expected 3)")
    pk = {}
    pp = os.path.join(a.kstar, "peak.tsv")
    if os.path.exists(pp):
        with open(pp) as fh:
            pk = dict(x.rstrip("\n").split("\t", 1) for x in fh if "\t" in x)
    line("peak_method", pk.get("method") == "genomescope", f"peak {pk.get('haploid_peak', 'NA')} by {pk.get('method', 'NA')} "
         "(expected genomescope on real data)", warn=True)
    for r in g.itertuples(index=False):
        rec = pd.read_csv(os.path.join(a.kstar, f"records_{r.uset[4:]}.tsv.gz"), sep="\t", dtype={"key": str})
        cok, cdet = calibration(rec, a.calib_min)
        line(f"calibration_{r.uset}", cok, cdet)
        sh = r.share_giab_tp_merfin_flagged_kstar
        raw = gs[(gs["uset"] == r.uset) & (gs["metric"] == "either_raw") & (gs["threshold"] == 1.0)]["share_tp_flagged"]
        fp = r.share_giab_fp_merfin_flagged_kstar
        line(f"tp_flagged_{r.uset}", bool(sh == sh and sh <= a.tp_max),
             f"Merfin-supported GIAB TP flagged by kstar {r.giab_tp_merfin_flagged_kstar} = {sh:.4f} (bound {a.tp_max}); "
             f"unrounded {float(raw.iloc[0]) if len(raw) else float('nan'):.4f}; supported FP flagged {fp:.4f}; "
             f"supported TP / FP without k* {r.giab_tp_merfin_kstar_na} / {r.giab_fp_merfin_kstar_na}", warn=cok)
    print(f"warnings\t{nwarn}\tfailed\t{nfail}")
    print("KSTAR_VERIFY_OK" if nfail == 0 else "KSTAR_VERIFY_FAILED")
    sys.exit(1 if nfail else 0)


def mode_state(a):
    samples = pd.read_csv(os.path.join(a.sheet_dir, "samples.tsv"), sep="\t", dtype=str, keep_default_na=False)
    me = samples[(samples["platform"] == a.platform) & (samples["sample"] == a.sample)]
    if len(me) != 1:
        print("not_in_sheet")
        return
    print(sample_states(me, a.results, a.state)[(a.platform, a.sample)])


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="mode", required=True)
    p = sp.add_parser("peak")
    p.add_argument("--hist", required=True)
    p.add_argument("--genomescope", default="")
    p.add_argument("--out", required=True)
    p = sp.add_parser("windows")
    p.add_argument("--keys", required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--k", type=int, default=31)
    p.add_argument("--fasta", required=True)
    p.add_argument("--index", required=True)
    for name in ("kmers", "compare"):
        p = sp.add_parser(name)
        p.add_argument("--reads-wig", required=True)
        p.add_argument("--ref-wig", required=True)
        p.add_argument("--peak", required=True)
        p.add_argument("--chunk", type=int, default=5000000)
        if name == "kmers":
            p.add_argument("--index", required=True)
            p.add_argument("--spot-fasta", required=True)
            p.add_argument("--spot-dump", required=True)
            p.add_argument("--out", required=True)
        else:
            p.add_argument("--dump", required=True)
            p.add_argument("--partial", action="store_true")
    for name in ("tables", "state"):
        p = sp.add_parser(name)
        p.add_argument("--sheet-dir", required=True)
        p.add_argument("--results", required=True)
        p.add_argument("--state", default="")
        p.add_argument("--platform", required=True)
        p.add_argument("--sample", default="HG002")
        if name == "tables":
            p.add_argument("--sets", required=True)
            p.add_argument("--uv-out", required=True)
            p.add_argument("--kstats", required=True)
            p.add_argument("--out", required=True)
            p.add_argument("--kstar-min", type=float, default=1.0)
            p.add_argument("--haploid-xy", action="store_true")
            p.add_argument("--min-kmers", type=int, default=5)
    p = sp.add_parser("verify")
    p.add_argument("--kstar", required=True)
    p.add_argument("--hg002", required=True)
    p.add_argument("--platform", required=True)
    p.add_argument("--tp-max", type=float, default=0.15)
    p.add_argument("--calib-min", type=int, default=1000)
    a = ap.parse_args()
    {"peak": mode_peak, "windows": mode_windows, "kmers": mode_kmers, "compare": mode_compare, "tables": mode_tables,
     "verify": mode_verify, "state": mode_state}[a.mode](a)


if __name__ == "__main__":
    main()
