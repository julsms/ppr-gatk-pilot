#!/bin/bash
# Graph-space alignment metrics of the PPR call set: vg stats -a on the GAM of the PPR pipeline, and a copy of its giraffe log,
# in the folder that mapping_summary.py and the figure scripts read.
# Usage: bash ppr_mapping_stats.sh <PPR GAM> <PPR giraffe log (_align.err)> <folder, for example Illumina_PPR_mapping_stats>
set -euo pipefail

GAM=$1
ERR=$2
M=$3/mapping
THREADS=${THREADS:-4}
mkdir -p $M

vg stats -p $THREADS -a $GAM > $M/vg_stats_a.txt
cp $ERR $M/
