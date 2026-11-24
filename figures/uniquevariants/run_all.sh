#!/bin/bash
# Figure 4 and Supplementary Figures S4 and S13-S18 from the outputs of analyses/uniquevariants: the data tables first, then the figures.
# Usage: UV_DIR=<output folder of analyses/uniquevariants> bash run_all.sh
set -euo pipefail

UV_DIR=$(cd "${UV_DIR:-uniquevariants}" && pwd)
export UV_DIR
export RELATED=${RELATED:-6}        # one of the two first-degree relatives, left out of the singletons
cd "$(dirname "$0")"

python3 uniquevariants_table.py
python3 figure4_panels.py
python3 figure4_combine.py
python3 supp_uniquevariants_why_gatk_missed.py
python3 supp_uniquevariants_hg002_precision.py
python3 supp_uniquevariants_frequency_classes.py
python3 supp_uniquevariants_1000g.py
python3 supp_uniquevariants_consequence_classes.py
python3 supp_uniquevariants_titv_hethom.py
python3 supp_uniquevariants_singletons.py
python3 titv_by_class.py "$UV_DIR/tables/r27_titv_hethom.tsv" titv_by_class.tsv
