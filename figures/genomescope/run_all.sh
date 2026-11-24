#!/bin/bash
# GenomeScope tables, the platform test of the call-set heterozygosity and the three GenomeScope supplementary figures.
# Run from the folder that holds genomescope/ (the outputs of analyses/genomescope) and uniquevariants/, or set GENOMESCOPE_DIR
# and PAPER_AGAIN_DIR.
# Usage: bash run_all.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)

python3 $HERE/genomescope_stats.py
python3 $HERE/het_from_calls.py
python3 $HERE/het_platform_test.py
python3 $HERE/supp_kmer_spectra.py
python3 $HERE/supp_genomescope_estimates.py
python3 $HERE/supp_genomescope_het_vs_length.py
