#!/bin/bash
# bench_chr20.sh (chromosomes 20-22) for the six call sets of Figure 1 (GATK, BWA-MEM + DeepVariant and PPR; Illumina and Aviti), one
# after the other, with nice; a call set with <name>/chr20_22/CHR20_22_DONE is skipped. Run it from the benchmark folder.
# STRATA=1 runs bench_chr20.sh with the strata (outputs <name>/chr20_22_strata/, marker CHR20_22_STRATA_DONE) and adds the three other
# graph set-ups of Figure 2 on both platforms (GRCh38-only, HPRC default, HPRC d9): 12 call sets.
# Usage: bash run_chr20_benchmarks.sh <list>      (<list>: one line per call set, the name and the query VCF separated by a tab)
#        STRATA=1 bash run_chr20_benchmarks.sh <list>
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
LIST=$(readlink -f $1)
export STRATA=${STRATA:-0}
TAG=chr20_22
NAMES=(Illumina_GATK_HaplotypeCaller_joint_vqsrbp_snp99.9_indel99.9 Illumina_GATK_Deepvariant_markdup Illumina_PPR_Deepvariant
       Aviti_GATK_HaplotypeCaller_joint_vqsrbp_snp99.9_indel99.9 Aviti_GATK_Deepvariant_markdup Aviti_PPR_Deepvariant)
if [[ $STRATA == 1 ]]; then
    TAG+=_strata
    for P in Illumina Aviti; do
        NAMES+=(${P}_giraffe_GRCh38only_Deepvariant ${P}_giraffe_HPRC_default_Deepvariant ${P}_giraffe_HPRC_d9_Deepvariant)
    done
fi
MARK=${TAG^^}
fail=0
for N in "${NAMES[@]}"; do
    [[ -f $N/$TAG/${MARK}_DONE ]] && continue
    VCF=$(awk -F'\t' -v n=$N '$1 == n {print $2}' $LIST)
    [[ -n $VCF ]] || { echo "no VCF for $N in $LIST"; fail=1; continue; }
    nice -n 10 bash $HERE/bench_chr20.sh $N $VCF || { echo "ERROR $N"; fail=1; }
done
n=0; for N in "${NAMES[@]}"; do [[ -f $N/$TAG/${MARK}_DONE ]] && n=$((n + 1)); done
if (( fail || n != ${#NAMES[@]} )); then echo ${MARK}_ALL_FAILED; exit 1; fi
echo ${MARK}_ALL_DONE
