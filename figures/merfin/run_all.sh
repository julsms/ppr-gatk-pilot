#!/bin/bash
# Data and figures of the k-mer support: merfin_table.py, then Figure 5 and the five Merfin supplementary figures.
# Run from the folder that holds kmers_validation/revision (the outputs of analyses/kmers), or set KV_DIR to it.
# Usage: bash run_all.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)

python3 $HERE/merfin_table.py
python3 $HERE/figure5_combine.py
python3 $HERE/supp_merfin_per_individual.py
python3 $HERE/supp_merfin_regions_cohort.py
python3 $HERE/supp_merfin_counts_by_region.py
python3 $HERE/supp_merfin_hg002_giab_classes.py
python3 $HERE/supp_merfin_hg002_regions_aviti.py
