#!/bin/bash
# Depth correlations and their figures: coverage_table.py, then Figure 6 and the depth-by-region supplementary figure.
# Run from the folder that holds qc/, concordance/ and kmers_validation/revision/ (analysis outputs), or set PAPER_AGAIN_DIR.
# Usage: bash run_all.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)

python3 $HERE/coverage_table.py
python3 $HERE/figure6.py
python3 $HERE/supp_depth_by_region.py
python3 $HERE/holm_by_variant_type.py
