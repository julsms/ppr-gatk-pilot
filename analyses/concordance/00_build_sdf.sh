#!/bin/bash
# The RTG SDF of the chr-named GRCh38 reference, the template of vcfeval (03_vcfeval_pair.sh).
# Usage: bash 00_build_sdf.sh <GRCh38 FASTA with chr names> ref/GRCh38_chr.sdf
set -euo pipefail

export RTG_MEM=${RTG_MEM:-12g}
mkdir -p $(dirname $2)
rtg format -o $2 $1
