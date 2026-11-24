#!/bin/bash
# One sample from a joint VCF (HG002 or HG001 for the benchmark): its column only, only the ALT alleles it carries,
# only the records where it carries one. Nothing is split or normalised and FILTER is kept.
# Usage: bash extract_sample.sh <joint VCF> <sample name> <output VCF>
set -euo pipefail

bcftools view -s $2 -a -c 1 -Oz -o $3 $1
bcftools index -f -t $3
