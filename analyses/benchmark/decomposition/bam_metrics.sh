#!/bin/bash
# Alignment metrics of one BAM, the same definitions for every set-up: primary records, mapped, properly paired, MAPQ >= 60,
# exact (NM:i:0 and a CIGAR of matches only) and gapless (no I or D, clipping allowed). Secondary and supplementary records are left out.
# KEEP_CONTIGS=chr1_22XY counts a mapped record only on chromosomes 1-22, X and Y (with or without the chr prefix); a record on another
# contig (chrM, unplaced) is counted as primary, not as mapped, like the giraffe BAMs, whose surjection keeps only these 24 paths. The
# default (empty) counts every contig. The run split the SAM lines over several copies of the same awk counter and summed the counts
# (AWK_WORKERS, a speed setting only); one awk gives the same counts.
# Usage: [KEEP_CONTIGS=chr1_22XY] bash bam_metrics.sh <arm> <BAM> <output folder>
#        (writes <output folder>/<arm>.tsv and <output folder>/BAM_METRICS_<arm>_DONE)
set -euo pipefail

ARM=$1
BAM=$2
D=$3
mkdir -p $D
case "${KEEP_CONTIGS:-}" in
    '') KEEP_RE='' ;;
    chr1_22XY) KEEP_RE='^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$' ;;
    *) echo "bad KEEP_CONTIGS '$KEEP_CONTIGS' (empty or chr1_22XY)"; exit 1 ;;
esac

samtools view -@ 4 -F 0x900 $BAM | awk -F'\t' -v keep="$KEEP_RE" '
    { n++; f = $2 + 0
      if (int(f / 4) % 2 == 1) next
      if (keep != "" && $3 !~ keep) next
      m++
      if (int(f / 2) % 2 == 1) pp++
      if ($5 + 0 >= 60) q60++
      c = $6
      if (c !~ /[ID]/) gl++
      nm = -1
      for (i = 12; i <= NF; i++) if (substr($i, 1, 5) == "NM:i:") { nm = substr($i, 6) + 0; break }
      if (nm < 0) nonm++
      else if (nm == 0 && c ~ /^[0-9]+M$/) ex++ }
    END { printf "primary\t%d\nmapped\t%d\nproperly_paired\t%d\nmapq60\t%d\nexact\t%d\ngapless\t%d\nno_NM_tag\t%d\n", n, m, pp, q60, ex, gl, nonm }' \
    > $D/$ARM.tsv
touch $D/BAM_METRICS_${ARM}_DONE
