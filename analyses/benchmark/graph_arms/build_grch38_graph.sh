#!/bin/bash
# The GRCh38-only graph: the 179 GRCh38 paths of the HPRC v1.1 default graph as FASTA, indexed for giraffe with vg autoindex.
# The paths keep their names (GRCh38#0#chr1 ...), so surject and DeepVariant run as for the HPRC graphs.
# Usage: bash build_grch38_graph.sh <hprc-v1.1-mc-grch38.gbz> <output folder>     (writes grch38.giraffe.gbz, grch38.dist, grch38.min)
set -euo pipefail

GBZ=$(readlink -f $1)
OUT=$2
THREADS=${THREADS:-12}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p $OUT/tmp
cd $OUT
MOUNTS=$(for d in $PWD $(dirname $GBZ); do echo "-v $d:$d"; done | sort -u)
VG="docker run --rm -u $(id -u):$(id -g) -e HOME=/tmp $MOUNTS -w $PWD quay.io/vgteam/vg:v1.52.0@sha256:795887f8381623d0242743b5fbb813e0fcc1b0fdd8dae04f335e0798f1e454e5"

# the FASTA must be the one of the run (grch38_graph.md5 also lists the index files of the run, for reference only)
$VG vg paths -x $GBZ -F -S GRCh38 > grch38_from_hprc.fa
grep grch38_from_hprc.fa $HERE/grch38_graph.md5 | md5sum -c -
$VG vg autoindex --workflow giraffe -r grch38_from_hprc.fa -p grch38 -t $THREADS -T tmp > autoindex.log 2>&1
rm -rf tmp
