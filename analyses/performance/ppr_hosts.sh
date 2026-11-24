#!/bin/bash
# PPR against GATK on the same individuals (ppr_v3_matched_gatk_summary.tsv), and the PPR time, CPU and memory by host and
# threads per step (ppr_v3_hosts.tsv, ppr_v3_hosts_summary.tsv), from the three tables of the trace scripts.
# Usage: bash ppr_hosts.sh [folder with ppr_v3_resources_per_sample.tsv and the three .csv tables]
set -uo pipefail
B=${1:-.}
cd "$B"

awk -v OFS='\t' '
    function add(l, v) { return (index("," l ",", "," v ",") ? l : (l == "" ? v : l "," v)) }
    FILENAME == "ppr_v3_resources_per_sample.tsv" {
        if (FNR == 1) { for (i = 1; i <= NF; i++) c[$i] = i; next }
        k = $c["platform"] "|" $c["sample"]
        host[k] = $c["host"]; prof[k] = $c["profile"]; heavy[k] = $c["cpus_heavy_step"]; dv[k] = $c["cpus_deepvariant"]
        img[k] = ($c["deepvariant_container"] ~ /-gpu$/) ? "GPU" : (($c["deepvariant_container"] ~ /deepvariant/) ? "CPU" : "NA")
        ord[++n] = k; next
    }
    FNR == 1 { next }
    FILENAME == "total_time_results.csv"   { v["time"   SUBSEP $1 "|" $2] = $3; next }
    FILENAME == "cpu_minutes_results.csv"  { v["cpu"    SUBSEP $1 "|" $2] = $3; next }
    FILENAME == "memory_usage_results.csv" { v["maxrss" SUBSEP $1 "|" $2] = $3; next }
    END {
        H = "ppr_v3_hosts.tsv"; S = "ppr_v3_hosts_summary.tsv"; M = "ppr_v3_matched_gatk_summary.tsv"
        print "platform", "sample", "host", "profile", "cpus_heavy_step", "cpus_deepvariant", "deepvariant_image", "total_time_h", "cpu_minutes", "max_peak_rss_gb" > H
        nm = split("time cpu maxrss", metrics, " ")
        for (i = 1; i <= n; i++) {
            k = ord[i]; split(k, a, "|")
            p = (a[1] == "illumina") ? "PPR_Illumina" : "PPR_Aviti"; g = (a[1] == "illumina") ? "GATK_Illumina" : "GATK_Aviti"
            # the GATK Aviti trace label of each PPR sample (named by individual)
            gs = a[2]; if (a[1] == "aviti" && gs == "12") gs = "21"; else if (a[1] == "aviti" && gs == "21") gs = "12"
            pk = p "|" a[2]; gk = g "|" gs
            print a[1], a[2], host[k], prof[k], heavy[k], dv[k], img[k], v["time" SUBSEP pk], v["cpu" SUBSEP pk], v["maxrss" SUBSEP pk] > H
            h = host[k]; hn[h]++; hp[h] = add(hp[h], prof[k]); hh[h] = add(hh[h], heavy[k]); hd[h] = add(hd[h], dv[k])
            hi[h] = add(hi[h], img[k])
            ht[h] += v["time" SUBSEP pk]; hc[h] += v["cpu" SUBSEP pk]
            if (v["maxrss" SUBSEP pk] + 0 > hm[h]) hm[h] = v["maxrss" SUBSEP pk] + 0
            for (j = 1; j <= nm; j++) {
                m = metrics[j]
                if ((m SUBSEP pk) in v && (m SUBSEP gk) in v) {
                    mk = a[1] SUBSEP m; np[mk]++; sp[mk] += v[m SUBSEP pk]; sg[mk] += v[m SUBSEP gk]
                }
            }
        }
        print "host", "profile", "n_samples", "cpus_heavy_step", "cpus_deepvariant", "deepvariant_image", "mean_total_time_h", "mean_cpu_minutes", "max_peak_rss_gb" > S
        close(S)
        for (h in hn) printf "%s\t%s\t%d\t%s\t%s\t%s\t%.2f\t%.0f\t%.1f\n", h, hp[h], hn[h], hh[h], hd[h], hi[h], ht[h] / hn[h], hc[h] / hn[h], hm[h] | "sort >> " S
        close("sort >> " S)
        print "platform", "metric", "n_pairs", "mean_PPR_v3", "mean_GATK", "ratio_PPR_over_GATK" > M
        split("illumina aviti", pls, " ")
        for (q = 1; q <= 2; q++) for (j = 1; j <= nm; j++) {
            mk = pls[q] SUBSEP metrics[j]
            if (np[mk] > 0) printf "%s\t%s\t%d\t%.2f\t%.2f\t%.3f\n", pls[q], metrics[j], np[mk], sp[mk] / np[mk], sg[mk] / np[mk], sp[mk] / sg[mk] > M
        }
    }
' FS='\t' ppr_v3_resources_per_sample.tsv FS=',' total_time_results.csv cpu_minutes_results.csv memory_usage_results.csv

echo "== per host (ppr_v3_hosts_summary.tsv)"; column -t -s $'\t' ppr_v3_hosts_summary.tsv
echo "== PPR vs GATK on the same samples (ppr_v3_matched_gatk_summary.tsv)"; column -t -s $'\t' ppr_v3_matched_gatk_summary.tsv
echo "included samples: $(($(wc -l < ppr_v3_hosts.tsv) - 1))"
