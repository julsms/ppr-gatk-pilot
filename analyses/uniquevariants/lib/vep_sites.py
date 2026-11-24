#!/usr/bin/env python3
"""vep_sites.py - VEP on the distinct alleles of all unique sets. sites: the sites-only VCF of the alleles (keys tables with a
.DONE file); parse: one row per allele, the fields read by name from the CSQ header (transcript: MANE Select, else canonical,
else the first); merge: one row per allele over the finished batches. novel = not in dbSNP b157 (custom, exact allele), gnomAD
genomes or 1000 Genomes. An allele is in dbSNP when the custom annotation gives the ID of the matched record (dbSNP) or its RS
(dbSNP_RS): INFO/RS is a 32-bit integer in htslib, so the prepared dbSNP file has no RS for rs numbers above 2,147,483,647.
"""
import argparse
import glob
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uvlib  # noqa: E402

FIELDS = ["Consequence", "IMPACT", "SYMBOL", "Gene", "BIOTYPE", "MANE_SELECT", "CANONICAL", "VARIANT_CLASS",
          "Existing_variation", "dbSNP", "dbSNP_RS", "AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF", "gnomADe_AF", "gnomADg_AF",
          "MAX_AF", "CLIN_SIG", "ClinVar_CLNSIG"]
KG = ["AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF"]


def any_pos(vals):
    for v in vals:
        for x in str(v).split("&"):
            try:
                if float(x) > 0:
                    return True
            except ValueError:
                pass
    return False


def finished_batches(bdir):
    """the batch tables of the finished batches (BATCH_DONE), oldest first (the folder names are date stamps)"""
    if not bdir:
        return []
    return [f for f in sorted(glob.glob(os.path.join(bdir, "*", "vep_sites.batch.tsv.gz")))
            if os.path.exists(os.path.join(os.path.dirname(f), "BATCH_DONE"))]


def sites(a):
    keys = set()
    for root in a.roots.split(","):
        for f in sorted(glob.glob(os.path.join(root, "*", "*.keys.tsv.gz"))):
            if not os.path.exists(f[: -len(".keys.tsv.gz")] + ".DONE"):
                continue
            k = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False, usecols=["CHROM", "POS", "REF", "ALT"])
            keys |= set(uvlib.make_key(k))
    order = {}
    hdr = ["##fileformat=VCFv4.2"]
    for i, ln in enumerate(open(a.fai)):
        c, n = ln.split("\t")[:2]
        order[c] = i
        hdr.append(f"##contig=<ID={c},length={n}>")
    done = set()
    for f in finished_batches(a.batches):
        done |= set(pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False, usecols=["key"])["key"])
    new = keys - done
    uvlib.write_tsv(pd.DataFrame({"key": sorted(keys)}), a.out_prefix + ".all_keys.tsv.gz")

    def sk(k):
        c, p, r, x = k.split(":")
        return (order[c], int(p), r, x)
    with open(a.out_prefix + ".sites.vcf", "w") as f:
        f.write("\n".join(hdr) + "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for k in sorted(new, key=sk):
            c, p, r, x = k.split(":")
            f.write(f"{c}\t{p}\t.\t{r}\t{x}\t.\t.\t.\n")
    print(f"keys {len(keys)}")
    print(f"new {len(new)}")


def parse(a):
    fmt, vep_line, rows = None, ".", []
    for line in uvlib.bcftools_lines(a.bcftools, ["view", a.vep_vcf]):
        if line.startswith("##INFO=<ID=CSQ"):
            fmt = line.split("Format: ")[1].split('"')[0].split("|")
            continue
        if line.startswith("##VEP="):
            vep_line = line.strip()
            continue
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        csq = ""
        for kv in f[7].split(";"):
            if kv.startswith("CSQ="):
                csq = kv[4:]
        rows.append((f"{f[0]}:{f[1]}:{f[3]}:{f[4]}", csq))
    if fmt is None:
        sys.exit("vep_sites.py parse: no CSQ header in the VEP output")
    ix = {n: i for i, n in enumerate(fmt)}
    missing = [n for n in FIELDS if n not in ix]
    out = []
    for key, csq in rows:
        tr = [t.split("|") for t in csq.split(",")] if csq else []
        pick = None
        for t in tr:
            if "MANE_SELECT" in ix and t[ix["MANE_SELECT"]]:
                pick = t
                break
        if pick is None:
            for t in tr:
                if "CANONICAL" in ix and t[ix["CANONICAL"]] == "YES":
                    pick = t
                    break
        if pick is None and tr:
            pick = tr[0]
        rec = {"key": key}
        for n in FIELDS:
            rec[n] = pick[ix[n]] if (pick is not None and n in ix) else ""
        rec["consequence_class"] = rec["Consequence"].split("&")[0] if rec["Consequence"] else "none"
        rec["in_dbSNP"] = int(rec["dbSNP"] not in ("", ".") or rec["dbSNP_RS"] not in ("", "."))
        rec["in_existing_rs"] = int("rs" in rec["Existing_variation"])
        rec["in_gnomADg"] = int(any_pos([rec["gnomADg_AF"]]))
        rec["in_gnomADe"] = int(any_pos([rec["gnomADe_AF"]]))
        rec["in_1000G"] = int(any_pos([rec[x] for x in KG]))
        rec["novel"] = int(not (rec["in_dbSNP"] or rec["in_gnomADg"] or rec["in_1000G"]))
        rec["novel_old"] = int(rec["novel"] and not rec["in_gnomADe"])
        out.append(rec)
    df = pd.DataFrame(out, columns=["key"] + FIELDS + ["consequence_class", "in_dbSNP", "in_existing_rs", "in_gnomADg",
                                                       "in_gnomADe", "in_1000G", "novel", "novel_old"])
    uvlib.write_tsv(df, a.out)
    with open(a.out + ".info.tsv", "w") as f:
        f.write(f"vep_header\t{vep_line}\nrecords\t{len(df)}\ncsq_fields\t{len(fmt)}\n")
        f.write(f"missing_fields\t{','.join(missing) or 'none'}\n")
        # the versions in the ##VEP line (VEP, API, cache, assembly, 1000genomes, dbSNP, gnomADe, gnomADg, ClinVar, ...)
        for k, v in re.findall(r'([A-Za-z0-9_.-]+)=("[^"]*"|\S+)', vep_line.replace("##VEP=", "VEP=", 1)):
            f.write(f"vep_db_{k}\t{v.strip(chr(34))}\n")
    print(f"parse: {len(df)} records; missing CSQ fields: {','.join(missing) or 'none'}")


def merge(a):
    parts = [pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False) for f in finished_batches(a.batches)]
    allk = pd.read_csv(a.all_keys, sep="\t", dtype=str, keep_default_na=False)
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["key"])
    df = df[df["key"].isin(set(allk["key"]))]
    dup = int(df["key"].duplicated().sum())
    # the same key in two batches must carry the same annotation; the latest batch is kept
    conflict = int((df.drop_duplicates().groupby("key").size() > 1).sum()) if dup else 0
    df = df.drop_duplicates("key", keep="last")
    miss = len(set(allk["key"]) - set(df["key"]))
    with open(a.out + ".checks.tsv", "w") as f:
        f.write(f"check\tvalue\ndistinct_keys\t{len(allk)}\nvep_rows\t{len(df)}\nkeys_in_two_batches\t{dup}\n"
                f"keys_with_differing_rows\t{conflict}\nkeys_without_row\t{miss}\n")
    if conflict or miss or len(df) != len(allk):
        sys.exit(f"vep_sites.py merge: {len(df)} rows for {len(allk)} keys ({conflict} keys with differing rows in two batches, "
                 f"{miss} without a row)")
    uvlib.write_tsv(df.sort_values("key"), a.out)
    print(f"merge: {len(df)} rows = {len(allk)} distinct keys")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["sites", "parse", "merge"])
    ap.add_argument("--roots")
    ap.add_argument("--fai")
    ap.add_argument("--out-prefix")
    ap.add_argument("--vep-vcf")
    ap.add_argument("--out")
    ap.add_argument("--batches")
    ap.add_argument("--all-keys")
    ap.add_argument("--bcftools", default=os.environ.get("BCFTOOLS", "bcftools"))
    a = ap.parse_args()
    {"sites": sites, "parse": parse, "merge": merge}[a.mode](a)


if __name__ == "__main__":
    main()
