#!/bin/bash
# Supplementary Figure S28 and the values of Supplementary Note 2: performance_table.py, then supp_performance.py.
# Set TIMELINE_DIR to the folder with the tables of analyses/performance.
# Usage: TIMELINE_DIR=<folder> bash run_all.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)

python3 $HERE/performance_table.py
python3 $HERE/supp_performance.py
