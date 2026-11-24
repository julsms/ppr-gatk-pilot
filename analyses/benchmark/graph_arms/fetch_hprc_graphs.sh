#!/bin/bash
# Downloads the HPRC v1.1 GRCh38 Minigraph-Cactus graphs with their giraffe indexes from the public HPRC bucket: the default graph
# (hprc-v1.1-mc-grch38, all haplotypes, the graph PPR samples from) and the d9 graph (nodes on at least 9 of the 90 haplotypes).
# Usage: bash fetch_hprc_graphs.sh <output folder>
set -euo pipefail

OUT=$1
HERE=$(cd "$(dirname "$0")" && pwd)
S3=s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38
mkdir -p $OUT
cd $OUT

for f in hprc-v1.1-mc-grch38.gbz hprc-v1.1-mc-grch38.dist hprc-v1.1-mc-grch38.min \
         hprc-v1.1-mc-grch38.d9.gbz hprc-v1.1-mc-grch38.d9.dist hprc-v1.1-mc-grch38.d9.min; do
    aws s3 cp --no-sign-request --region us-west-2 --endpoint-url https://s3.us-west-2.amazonaws.com --only-show-errors $S3/$f $f
done
md5sum -c $HERE/hprc_graphs.md5
