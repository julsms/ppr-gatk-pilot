#!/bin/bash
# Cohort tables of the k-mer support (validation_tables.py) and Merfin against GIAB on HG002 (hg002_merfin_giab.py).
# Run after merfin_sample.sh on every sample of the sheet.
# Usage: bash tables.sh <output folder> <unique-variant output folder>
set -euo pipefail

OUT=$1
UV=$2                       # analyses/uniquevariants outputs: vep/vep_sites.tsv.gz and hg002/truth/
HERE=$(cd "$(dirname "$0")" && pwd)

python3 $HERE/validation_tables.py --sheet-dir $OUT/sheet --results $OUT/results --state $OUT/sample_state.tsv \
    --out $OUT/tables --vep $UV/vep/vep_sites.tsv.gz
python3 $HERE/hg002_merfin_giab.py --sheet-dir $OUT/sheet --results $OUT/results --state $OUT/sample_state.tsv \
    --sets $HERE/config/merfin_sets.tsv --uv-out $UV --out $OUT/hg002
