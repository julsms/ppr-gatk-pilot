# Writes every genotype unphased with the alleles sorted: 0|1, 1|0 and 1/0 become 0/1, 1|1 becomes 1/1, 1|. and ./1 become ./1.
# Only GT (the first FORMAT field) changes. A record without GT first stops with exit 1.
# Usage: bcftools view -Ov in.vcf.gz | awk -f lib/unphase_sort_gt.awk | bgzip > out.vcf.gz
BEGIN { FS = OFS = "\t" }
/^#/ { print; next }
{
    nrec++
    if ($9 !~ /^GT(:|$)/) { print "unphase_sort_gt.awk: GT is not the first FORMAT field in record " nrec > "/dev/stderr"; exit 1 }
    for (i = 10; i <= NF; i++) {
        c = index($i, ":")
        if (c > 0) { gt = substr($i, 1, c - 1); rest = substr($i, c) } else { gt = $i; rest = "" }
        n = split(gt, a, /[\/|]/)
        # insertion sort; "." sorts before any allele index
        for (j = 2; j <= n; j++) {
            v = a[j]; k = j - 1
            while (k >= 1 && key(a[k]) > key(v)) { a[k + 1] = a[k]; k-- }
            a[k + 1] = v
        }
        out = a[1]
        for (j = 2; j <= n; j++) out = out "/" a[j]
        $i = out rest
    }
    print
}
function key(x) { return (x == ".") ? -1 : x + 0 }
