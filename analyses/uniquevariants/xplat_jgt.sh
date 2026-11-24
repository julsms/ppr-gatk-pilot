#!/bin/bash
# Genotype concordance between Illumina and Aviti of four unique sets, computed as for the full call sets: PPR vs GATK, GATK vs PPR,
# PPR vs DV and DV vs PPR. For each individual, the unique set on Illumina (side A) is compared with the same set on Aviti (side B)
# with RTG vcfeval (the settings of 03_vcfeval_pair.sh), and 04_classify_count.py gives the classes and J_gt per region and type.
# The key tables lose their BK_lm column first, which 04 would read as a region. Outputs in uniquevariants/xplat_jgt/: prep/,
# vcfeval/, classify/ (pairs of config/xplat_pairs.tsv). A finished step (.DONE) is skipped on a new start.
# Usage: bash xplat_jgt.sh <sample>      (after extract_unique.sh on both platforms; SETS to run fewer sets)
set -euo pipefail

S=$1                                # the individual, as in the concordance sample map
CONC_OUT=${CONC_OUT:-concordance}   # outputs of the concordance step: ref/
UV_OUT=${UV_OUT:-uniquevariants}
SDF=${SDF:-$CONC_OUT/ref/GRCh38_chr.sdf}
THREADS=${THREADS:-4}
export RTG_MEM=${RTG_MEM:-12g}
HERE=$(cd "$(dirname "$0")" && pwd)
CONC=$(cd "$HERE/../concordance" && pwd)
export PAIRS=$HERE/config/xplat_pairs.tsv
O=$UV_OUT/xplat_jgt

for U in ${SETS:-PPRvsGATK GATKvsPPR PPRvsDV DVvsPPR}; do
    PAIR=${U}_illumina_aviti
    for P in illumina aviti; do
        x=$UV_OUT/unique/${U}_$P/$S
        [[ -f $x.DONE ]] || { echo "no unique set $x" >&2; exit 1; }
        k=$O/prep/${U}_$P/$S.keys.tsv.gz
        [[ -s $k ]] || python3 $HERE/lib/keys_without_bk_lm.py $x.keys.tsv.gz $k
    done
    V=$O/vcfeval/$PAIR/$S
    if [[ ! -f $V.DONE ]]; then
        mkdir -p $O/vcfeval/$PAIR
        rtg vcfeval -b $UV_OUT/unique/${U}_illumina/$S.vcf.gz -c $UV_OUT/unique/${U}_aviti/$S.vcf.gz -t $SDF -o $V \
            --output-mode=ga4gh --ref-overlap -T $THREADS
        touch $V.DONE                   # 04_classify_count.py reads this file as the sign that vcfeval has finished
    fi
    [[ -f $O/classify/$PAIR/$S.DONE ]] || python3 $CONC/04_classify_count.py --pair $PAIR --patient $S \
        --a-keys $O/prep/${U}_illumina/$S.keys.tsv.gz --b-keys $O/prep/${U}_aviti/$S.keys.tsv.gz --vcfeval $V \
        --out-prefix $O/classify/$PAIR/$S
done
