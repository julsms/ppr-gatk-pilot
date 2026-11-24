#!/usr/bin/env python3
"""truth_preclean.py - the GIAB truth without the records whose REF differs from the reference (bcftools norm --check-ref).
GIAB v5.0q masks the false duplications of GRCh38 with N, so one deletion (chr9, outside the benchmark bed) has N in its REF and
the preparation stops on it. Every other line is kept unchanged; above --max-in-bb excluded records inside the bed it stops.
Outputs: <out>/<name>.refclean.vcf.gz (+ .tbi), the excluded records, excluded.tsv, summary.tsv, PRECLEAN_DONE
"""
import argparse
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
import time

CODE = os.path.abspath(__file__)
NON_ACGTN = re.compile(r"[^ACGTN]")


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def stamp(path):
    st = os.stat(path)
    return f"{st.st_size} {int(st.st_mtime)}"


def read_kv(path):
    out = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t", 1)
            if len(p) == 2:
                out[p[0]] = p[1]
    return out


def write_kv(path, rows):
    with open(path, "w") as f:
        for k, v in rows:
            f.write(f"{k}\t{v}\n")


def canon(seq):
    """the comparison of bcftools norm --check-ref: no case, IUPAC codes (any letter other than ACGTN) as N"""
    return NON_ACGTN.sub("N", seq.upper())


class Fasta:
    """reads reference bases through the .fai offsets (no whole-chromosome load)"""

    def __init__(self, path):
        self.f = open(path, "rb")
        self.idx = {}
        with open(path + ".fai") as fai:
            for line in fai:
                c, ln, off, lb, lw = line.split("\t")[:5]
                self.idx[c] = (int(ln), int(off), int(lb), int(lw))

    def fetch(self, chrom, pos, n):
        """n bases from 1-based pos (fewer at the contig end)"""
        ln, off, lb, lw = self.idx[chrom]
        s, e = pos - 1, min(pos - 1 + n, ln)
        if s >= e:
            return ""
        b = off + (s // lb) * lw + s % lb
        last = off + ((e - 1) // lb) * lw + (e - 1) % lb
        self.f.seek(b)
        return self.f.read(last - b + 1).decode().replace("\n", "").replace("\r", "")


def run(cmd, **kw):
    return subprocess.run(cmd, **kw)


def norm_lines(stderr_text):
    """the 'Lines total/split/joined/realigned/skipped' counts of bcftools norm"""
    for line in stderr_text.splitlines():
        if line.startswith("Lines"):
            return [int(x) for x in line.split("\t")[1].split("/")]
    return None


def var_type(ref, alts):
    if len(alts) > 1:
        return "MULTI"
    a = alts[0]
    if len(a) == len(ref):
        return "SNV" if len(ref) == 1 else "MNP"
    return "DEL" if len(a) < len(ref) else "INS"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--truth", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--bb", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--bgzip", default="bgzip")
    ap.add_argument("--tabix", default="tabix")
    ap.add_argument("--max-in-bb", type=int, default=10)
    a = ap.parse_args()
    t0 = time.time()

    for f in (a.truth, a.ref, a.ref + ".fai", a.bb):
        if not os.path.isfile(f) or os.path.getsize(f) == 0:
            sys.exit(f"truth_preclean: missing {f}")
    name = re.sub(r"\.vcf(\.gz)?$", "", os.path.basename(a.truth))
    out = os.path.abspath(a.out)
    done = os.path.join(out, "PRECLEAN_DONE")
    cleaned = os.path.join(out, f"{name}.refclean.vcf.gz")
    excluded_vcf = os.path.join(out, f"{name}.ref_mismatch_excluded.vcf.gz")
    bcf_version = run([a.bcftools, "--version"], stdout=subprocess.PIPE, text=True, check=True).stdout.splitlines()[0]

    src_md5 = md5_file(a.truth)
    fp_lines = [f"truth {os.path.abspath(a.truth)} {stamp(a.truth)} {src_md5}",
                f"ref {os.path.abspath(a.ref)} {stamp(a.ref)} {md5_file(a.ref + '.fai')}",
                f"bb {os.path.abspath(a.bb)} {md5_file(a.bb)}",
                f"max_in_bb {a.max_in_bb}", f"bcftools {bcf_version}", f"code {md5_file(CODE)}"]
    fp = hashlib.md5("\n".join(fp_lines).encode()).hexdigest()

    if os.path.exists(done):
        d = read_kv(done)
        if d.get("fingerprint") == fp and os.path.exists(cleaned) and md5_file(cleaned) == d.get("cleaned_md5"):
            print(f"truth_preclean: {done} is current: {d.get('records_excluded')} records excluded "
                  f"({d.get('excluded_in_bb')} inside BB); cleaned truth {cleaned}")
            print("PRECLEAN_CURRENT")
            return 0
        print(f"truth_preclean: {done} has another fingerprint (the source, the reference, the BB bed or the code changed), or "
              f"the cleaned file changed; nothing written. Move {out} aside and start again.")
        print("PRECLEAN_STALE")
        return 65

    tmp = os.path.join(out, "tmp")
    os.makedirs(tmp, exist_ok=True)

    # 1. the REF_MISMATCH warnings of bcftools norm --check-ref w
    err1 = os.path.join(tmp, "norm_check.stderr.txt")
    with open(err1, "w") as ef:
        r = run([a.bcftools, "norm", "-f", a.ref, "--check-ref", "w", "-Ou", "-o", os.devnull, a.truth], stderr=ef)
    if r.returncode != 0:
        sys.exit(f"truth_preclean: bcftools norm --check-ref w failed (exit {r.returncode}; {err1})")
    text = open(err1).read()
    lines1 = norm_lines(text)
    if lines1 is None:
        sys.exit(f"truth_preclean: no line count in {err1}")
    warn = {}
    n_warn = n_nonacgtn = 0
    for line in text.splitlines():
        p = line.split("\t")
        if p[0] == "REF_MISMATCH":
            n_warn += 1
            warn[(p[1], int(p[2]), p[3].upper())] = 0
        elif p[0] == "NON_ACGTN_REF":
            n_nonacgtn += 1
    if len(warn) != n_warn:
        sys.exit(f"truth_preclean: {n_warn} warnings name {len(warn)} distinct records")

    # 2. one pass over the truth: excluded records apart, every other line unchanged
    tag = (f'##truth_preclean=<Source="{os.path.abspath(a.truth)}",SourceMD5={src_md5},Reference="{os.path.abspath(a.ref)}",'
           f'Removed={n_warn},Description="records whose REF differs from the reference (bcftools norm --check-ref), removed '
           f'before the preparation of the truth">\n')
    tmp_clean = os.path.join(tmp, os.path.basename(cleaned))
    tmp_excl = os.path.join(tmp, os.path.basename(excluded_vcf))
    n_in = n_kept = 0
    excl = []
    with gzip.open(a.truth, "rt") as fin, open(tmp_clean, "wb") as fc, open(tmp_excl, "wb") as fx:
        pc = subprocess.Popen([a.bgzip, "-c"], stdin=subprocess.PIPE, stdout=fc, text=True)
        px = subprocess.Popen([a.bgzip, "-c"], stdin=subprocess.PIPE, stdout=fx, text=True)
        for line in fin:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    pc.stdin.write(tag)
                    px.stdin.write(tag)
                pc.stdin.write(line)
                px.stdin.write(line)
                continue
            n_in += 1
            p = line.split("\t", 10)
            k = (p[0], int(p[1]), p[3].upper())
            if k in warn:
                warn[k] += 1
                excl.append(p)
                px.stdin.write(line)
            else:
                n_kept += 1
                pc.stdin.write(line)
        pc.stdin.close()
        px.stdin.close()
        if pc.wait() != 0 or px.wait() != 0:
            sys.exit("truth_preclean: bgzip failed")

    # 3. checks
    bad = [f"{k[0]}:{k[1]} matched {v} records" for k, v in warn.items() if v != 1]
    if bad:
        sys.exit("truth_preclean: a warning does not match exactly one record: " + "; ".join(bad[:5]))
    if n_in != n_kept + len(excl) or n_in != lines1[0]:
        sys.exit(f"truth_preclean: records in {n_in}, kept {n_kept}, excluded {len(excl)}, bcftools norm lines {lines1[0]}")
    fa = Fasta(a.ref)
    rows = []
    for p in excl:
        chrom, pos, ref, alts = p[0], int(p[1]), p[3], p[4].split(",")
        seq = fa.fetch(chrom, pos, len(ref)) if chrom in fa.idx else ""
        if canon(seq) == canon(ref):
            sys.exit(f"truth_preclean: {chrom}:{pos}: bcftools warned, but the REF equals the reference when read here")
        nd = sum(x != y for x, y in zip(canon(seq), canon(ref))) if len(seq) == len(ref) else "."
        filt = p[6]
        gt = p[9].split(":")[0].rstrip("\n") if len(p) > 9 else "."
        rows.append({"CHROM": chrom, "POS": pos, "END": pos + len(ref) - 1, "REF_len": len(ref),
                     "ALT_len": max(len(x) for x in alts), "type": var_type(ref, alts), "FILTER": filt, "GT": gt,
                     "in_BB": "no", "ref_seq_N": canon(seq).count("N"), "truth_ref_N": canon(ref).count("N"),
                     "n_diff_bases": nd, "ref_seq_head": seq[:30].upper(), "truth_ref_head": ref[:30].upper()})
    run([a.tabix, "-f", "-p", "vcf", tmp_clean], check=True)
    err2 = os.path.join(tmp, "norm_check_cleaned.stderr.txt")
    with open(err2, "w") as ef:
        r = run([a.bcftools, "norm", "-f", a.ref, "--check-ref", "e", "-Ou", "-o", os.devnull, tmp_clean], stderr=ef)
    lines2 = norm_lines(open(err2).read())
    if r.returncode != 0 or lines2 is None or lines2[0] != n_kept:
        sys.exit(f"truth_preclean: bcftools norm --check-ref e on the cleaned truth: exit {r.returncode}, lines {lines2}, "
                 f"kept {n_kept} ({err2})")

    # 4. inside the benchmark bed
    in_bb = set()
    if excl:
        q = run([a.bcftools, "query", "-T", a.bb, "--targets-overlap", "variant", "-f", "%CHROM\t%POS\t%REF\t%ALT\n",
                 tmp_excl], stdout=subprocess.PIPE, text=True, check=True)
        for line in q.stdout.splitlines():
            c, ps, rf, al = line.split("\t")
            in_bb.add((c, int(ps), rf.upper()))
    for r_, p in zip(rows, excl):
        r_["in_BB"] = "yes" if (p[0], int(p[1]), p[3].upper()) in in_bb else "no"
    n_bb = sum(r_["in_BB"] == "yes" for r_ in rows)

    cols = ["CHROM", "POS", "END", "REF_len", "ALT_len", "type", "FILTER", "GT", "in_BB", "ref_seq_N", "truth_ref_N",
            "n_diff_bases", "ref_seq_head", "truth_ref_head"]
    with open(os.path.join(tmp, "excluded.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r_ in rows:
            f.write("\t".join(str(r_[c]) for c in cols) + "\n")
    byc = {}
    for r_ in rows:
        byc[(r_["CHROM"], r_["in_BB"])] = byc.get((r_["CHROM"], r_["in_BB"]), 0) + 1
    with open(os.path.join(tmp, "by_chrom.tsv"), "w") as f:
        f.write("CHROM\tin_BB\trecords\n")
        for (c, b), n in sorted(byc.items()):
            f.write(f"{c}\t{b}\t{n}\n")
    clean_md5 = md5_file(tmp_clean)
    summary = [("key", "value"), ("source", os.path.abspath(a.truth)), ("source_md5", src_md5), ("reference", a.ref),
               ("bb_bed", a.bb), ("records_input", n_in), ("records_kept", n_kept), ("records_excluded", len(excl)),
               ("excluded_in_bb", n_bb), ("excluded_outside_bb", len(excl) - n_bb), ("max_in_bb", a.max_in_bb),
               ("bcftools_ref_mismatch_warnings", n_warn), ("bcftools_non_acgtn_ref_warnings", n_nonacgtn),
               ("bcftools_norm_lines_input", lines1[0]), ("bcftools_norm_lines_cleaned", lines2[0]),
               ("cleaned_check_ref_e", "ok"), ("cleaned", cleaned), ("cleaned_md5", clean_md5),
               ("bcftools", bcf_version), ("code_md5", md5_file(CODE)), ("seconds", round(time.time() - t0, 1))]
    write_kv(os.path.join(tmp, "summary.tsv"), summary)
    msg = (f"truth_preclean: {n_in} truth records, {len(excl)} excluded (REF differs from the reference): {n_bb} inside BB, "
           f"{len(excl) - n_bb} outside; kept {n_kept}; by chromosome: "
           + (", ".join(f"{c} {b} {n}" for (c, b), n in sorted(byc.items())) or "none"))
    print(msg)
    if n_bb > a.max_in_bb:
        print(f"truth_preclean: {n_bb} excluded records inside BB, above --max-in-bb {a.max_in_bb}: the truth would change "
              f"inside the benchmark bed; no DONE written (tables in {tmp})")
        print("PRECLEAN_TOO_MANY_IN_BB")
        return 2

    # atomic publish
    for f in ("excluded.tsv", "by_chrom.tsv", "summary.tsv", "norm_check.stderr.txt", "norm_check_cleaned.stderr.txt"):
        shutil.move(os.path.join(tmp, f), os.path.join(out, f))
    shutil.move(tmp_excl, excluded_vcf)
    shutil.move(tmp_clean + ".tbi", cleaned + ".tbi")
    shutil.move(tmp_clean, cleaned)
    write_kv(done + ".tmp", [("date", time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())), ("fingerprint", fp),
                             ("source", os.path.abspath(a.truth)), ("source_md5", src_md5), ("cleaned", cleaned),
                             ("cleaned_md5", clean_md5), ("records_input", n_in), ("records_kept", n_kept),
                             ("records_excluded", len(excl)), ("excluded_in_bb", n_bb), ("summary", os.path.join(out, "summary.tsv"))])
    os.replace(done + ".tmp", done)
    print("PRECLEAN_DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
