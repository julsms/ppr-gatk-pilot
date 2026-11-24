# Renames the contigs of a Broad hg38 resource VCF (chr1, chrM, chrUn_GL000195v1, ...) to the Ensembl names of the reference
# (1, MT, GL000195.1, ...) and writes the contig lines of the reference .dict. Records on contigs that are not in the .dict are
# dropped and counted. The output is not sorted: pipe it to bcftools sort.
# Usage: gzip -cdf <in.vcf.gz> | awk -v DICT=<.dict> -v STATS=<stats.tsv> -v NOTE=<text> -f rename_contigs.awk | bcftools sort -Oz -o <out> -
BEGIN {
    FS = "\t"; OFS = "\t"
    n = 0
    while ((getline line < DICT) > 0) {
        if (line !~ /^@SQ/) continue
        sn = ""; ln = ""
        k = split(line, f, "\t")
        for (i = 1; i <= k; i++) {
            if (f[i] ~ /^SN:/) sn = substr(f[i], 4)
            if (f[i] ~ /^LN:/) ln = substr(f[i], 4)
        }
        n++; order[n] = sn; len[sn] = ln; indict[sn] = 1
    }
    close(DICT)
    if (n == 0) { print "no @SQ line in " DICT > "/dev/stderr"; exit 1 }
}
function ensembl_name(c,    parts, acc) {
    if (c ~ /_alt$/) return ""
    sub(/^chr/, "", c)
    if (c == "M") return "MT"
    if (index(c, "_") > 0) {          # chr1_KI270706v1_random -> KI270706.1; chrUn_GL000195v1 -> GL000195.1
        split(c, parts, "_"); acc = parts[2]; sub(/v/, ".", acc); return acc
    }
    return c
}
/^##contig=/ { next }
/^##/ { print; next }
/^#CHROM/ {
    for (i = 1; i <= n; i++) print "##contig=<ID=" order[i] ",length=" len[order[i]] ">"
    print "##ensembl_contig_renaming=" NOTE   # the key must not start with contig: htsjdk reads every ##contig... line as a contig line
    print
    next
}
{
    c = $1
    if (!(c in map)) { e = ensembl_name(c); map[c] = (e != "" && (e in indict)) ? e : "" }
    if (map[c] == "") { dropped[c]++; nd++; next }
    sub(/^[^\t]+/, map[c]); kept[c]++; nk++
    print
}
END {
    print "status\tcontig_in\tcontig_out\trecords" > STATS
    for (c in map) {
        if (map[c] != "") print "kept\t" c "\t" map[c] "\t" kept[c] > STATS
        else print "dropped\t" c "\t-\t" dropped[c] > STATS
    }
    print "total_kept\t-\t-\t" nk + 0 > STATS
    print "total_dropped\t-\t-\t" nd + 0 > STATS
}
