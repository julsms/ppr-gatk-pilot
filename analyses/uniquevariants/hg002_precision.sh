#!/bin/bash
# HG002 as an individual of the cohort: its three call sets prepared, the pairs GATK vs PPR, DV vs PPR and GATK vs DV (vcfeval
# and the classes), the five unique sets, and each full set against the prepared GIAB truth (whole analysis space).
# Usage: bash hg002_precision.sh <platform> <GATK joint VCF> <PPR HG002 VCF> <DV HG002 VCF>    (after hg002_truth.sh)
set -euo pipefail

P=$1
GATK=$2                             # the final joint VCF of the paper (25 samples, HG002 included)
PPR=$3                              # HG002 DeepVariant VCF of the PPR pipeline
DV=$4                               # HG002 DeepVariant VCF of BWA-MEM + DeepVariant
export REF=${REF:-GRCh38_chr.fa}    # GRCh38 with chr names (+ .fai), as in the concordance step
CONC_OUT=${CONC_OUT:-concordance}
SDF=${SDF:-$CONC_OUT/ref/GRCh38_chr.sdf}
THREADS=${THREADS:-8}
export RTG_MEM=${RTG_MEM:-12g}
H=${UV_OUT:-uniquevariants}/hg002
CONC=$(cd "$(dirname "$0")/../concordance" && pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
if [[ $P == illumina ]]; then GN=HG002_SRR14724532; PN=HG002Illumina_SRR14724532; DN=HG002_SRR14724532
else GN=HG002Aviti_SRR33551004; PN=HG002Aviti_SRR33551004; DN=HG002Aviti_SRR33551004; fi

# 1. HG002 from each call set, PASS only, prepared as the cohort
bash $HERE/prepare_callset.sh $DV $DN HG002 $H/prep/DVH002_$P/HG002 pass ensembl
bash $HERE/prepare_callset.sh $GATK $GN HG002 $H/prep/GATK_$P/HG002 pass ensembl
bash $HERE/prepare_callset.sh $PPR $PN HG002 $H/prep/PPRH002_$P/HG002 pass chr

# 2. the pairs, as in the concordance step (A baseline, B calls)
pair() {   # pair <name> <set A> <set B>
    mkdir -p $H/vcfeval/${1}_$P
    rtg vcfeval -b $H/prep/${2}_$P/HG002.vcf.gz -c $H/prep/${3}_$P/HG002.vcf.gz -t $SDF \
        -o $H/vcfeval/${1}_$P/HG002 --output-mode=ga4gh --ref-overlap -T $THREADS
    touch $H/vcfeval/${1}_$P/HG002.DONE      # 04_classify_count.py reads only finished runs
    OUT=$H PAIRS=$HERE/config/hg002_pairs.tsv python3 $CONC/04_classify_count.py --pair ${1}_$P --patient HG002
}
pair H_GATK_PPR GATK PPRH002
pair H_DV_PPR DVH002 PPRH002
pair H_GATK_DV GATK DVH002

# 3. the unique sets of HG002, as extract_unique.sh
extract() {   # extract <set> <side> <prepared set of that side> <pair>...
    local U=${1}_$P SIDE=$2 SET=$3 rec= cnt= pr
    shift 3
    for pr in "$@"; do
        rec=$rec,$H/classify/${pr}_$P/HG002.records.tsv.gz
        cnt=$cnt,$H/classify/${pr}_$P/HG002.counts.tsv
    done
    mkdir -p $H/unique/$U
    python3 $HERE/lib/extract_unique.py --records ${rec#,} --counts ${cnt#,} --side $SIDE \
        --keys $H/prep/${SET}_$P/HG002.keys.tsv.gz --vcf $H/prep/${SET}_$P/HG002.vcf.gz \
        --out-prefix $H/unique/$U/HG002 --uset $U --patient HG002
    touch $H/unique/$U/HG002.DONE       # vep.sh annotates only the sets with a .DONE file
}
extract PPRvsGATK B PPRH002 H_GATK_PPR
extract GATKvsPPR A GATK H_GATK_PPR
extract PPRvsDV B PPRH002 H_DV_PPR
extract DVvsPPR A DVH002 H_DV_PPR
extract PPRvsBOTH B PPRH002 H_GATK_PPR H_DV_PPR

# 4. each full set against the prepared truth; the benchmark regions (stratum BB) come from the two key tables
mkdir -p $H/truth
for SET in DVH002 GATK PPRH002; do
    rtg vcfeval -b $H/prep/GIAB_v5.0q/HG002.vcf.gz -c $H/prep/${SET}_$P/HG002.vcf.gz -t $SDF \
        -o $H/truth/${SET}_$P --output-mode=ga4gh --ref-overlap -T $THREADS
    python3 $HERE/lib/truth_status.py status --vcfeval $H/truth/${SET}_$P --keys $H/prep/${SET}_$P/HG002.keys.tsv.gz \
        --truth-keys $H/prep/GIAB_v5.0q/HG002.keys.tsv.gz --bb-stratum BB --out-prefix $H/truth/${SET}_$P
done
