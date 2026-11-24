#!/bin/bash
# The tables and figures of the concordance: Figure 3, Supplementary Figures S8-S12 and the genome and BB rows of Supplementary File 4.
# Usage: CONCORDANCE_DIR=<absolute path of the concordance output folder> bash run_all.sh
set -euo pipefail
cd "$(dirname "$0")"

# concordance_table.py checks its tests against those of the original run
mkdir -p $CONCORDANCE_DIR/stats
cp posthoc_long.tsv $CONCORDANCE_DIR/stats/
python3 concordance_table.py
python3 figure3_panels.py
python3 figure3_combine.py
python3 supp_concordance_crossplatform_tests.py
python3 supp_concordance_caller_vs_reference.py
python3 supp_concordance_all_pairs.py
python3 supp_concordance_across_platforms_lines.py
python3 supp_concordance_heatmap.py
python3 supp_table_concordance_counts.py
