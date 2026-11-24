#!/bin/bash
# Figure 7 from the hap.py summaries in data/HG002/ and data/HG001/ next to this script: figure7_table.py, then figure7.py.
# Usage: bash run_all.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)

python3 $HERE/figure7_table.py > $HERE/figure7_table.txt
python3 $HERE/figure7.py
