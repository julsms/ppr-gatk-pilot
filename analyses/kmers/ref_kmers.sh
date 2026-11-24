#!/bin/bash
# The GRCh38 k-mers (k = 31) for merfin -seqmers and for the k* analysis. Built once.
# Usage: bash ref_kmers.sh <reference fasta> <output .meryl>
set -euo pipefail

REF=$1
OUT=$2
THREADS=${THREADS:-16}
MEM=${MEM:-100}             # GB for meryl count

meryl count k=31 threads=$THREADS memory=$MEM $REF output $OUT
