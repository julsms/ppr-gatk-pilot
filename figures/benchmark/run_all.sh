#!/bin/bash
# The tables and figures of this folder (Figure 1, Supplementary Figures S1 to S3, Figure 2, Supplementary Figures S5 to S7 and
# Supplementary Table S1), in the order of the run; the outputs go next to the scripts and into figures/.
# Usage: bash run_all.sh <benchmark folder>
set -euo pipefail
export BENCHMARK_REVISION_DIR=$(readlink -f $1)
cd "$(dirname "$0")"

python3 benchmark_table.py
python3 figure1_panels.py
python3 figure1_combine.py
python3 supp_metrics_by_region.py
python3 supp_metrics_by_region.py --filter ALL
python3 chr20_table.py
python3 supp_chr20_pr.py
python3 chr20_compare.py
python3 graph_arms_table.py
python3 figure2_panels.py
python3 figure2_combine.py
python3 supp_graph_arms_f1.py
python3 supp_graph_arms_alignment.py
python3 supp_table_graph_setups_errors.py
