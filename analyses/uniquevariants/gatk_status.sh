#!/bin/bash
# Why GATK did not report the PPR-only alleles (vs GATK) of one individual: the raw GATK calls of the individual are added back in
# three steps (those removed by VQSR, then by the ExcessHet filter, then by any other filter), and each step is compared with the
# PPR calls. Run gatk_status_cohort.sh when every individual is done.
# Usage: bash gatk_status.sh <sample> <platform> <raw joint VCF> <final joint VCF>
set -euo pipefail

S=$1
P=$2
RAW=$3                              # merged_raw_variants.vcf.gz of the GATK joint calling (before ExcessHet and VQSR)
FINAL=$4                            # the final joint VCF of the paper (VQSR SNP 99.9 / INDEL 99.9), every FILTER value kept
CONC_OUT=${CONC_OUT:-concordance}
UV_OUT=${UV_OUT:-uniquevariants}
SDF=${SDF:-$CONC_OUT/ref/GRCh38_chr.sdf}     # the RTG SDF of the reference, as in the concordance step
THREADS=${THREADS:-8}
export RTG_MEM=${RTG_MEM:-12g}
CONC=$(cd "$(dirname "$0")/../concordance" && pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
R=$UV_OUT/reasons/new_$P/$S
U=$UV_OUT/prep
C=$CONC_OUT/classify/GATK_PPR_$P/$S
# the name of the individual in the joint VCFs (column <platform>_gatk of the concordance sample map)
G=$(awk -F'\t' -v s=$S -v c=${P}_gatk 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $h["patient"] == s {print $h[c]}' \
    $CONC/config/sample_map.tsv)

# 1. the raw calls and the final calls of the individual, every FILTER value kept, prepared as the call sets
bash $HERE/prepare_callset.sh $RAW $G $S $U/GATKRAW_$P/$S all ensembl
bash $HERE/prepare_callset.sh $FINAL $G $S $U/GATKALL_$P/$S all ensembl

# 2. the files that reason_states.py and reason_classes.py read
mkdir -p $R
printf '%s\t%s\n' pair_records $C.records.tsv.gz pair_counts $C.counts.tsv \
    ppr_vcf $CONC_OUT/prep/PPR_$P/$S.vcf.gz ppr_keys $CONC_OUT/prep/PPR_$P/$S.keys.tsv.gz \
    raw_vcf $U/GATKRAW_$P/$S.vcf.gz raw_keys $U/GATKRAW_$P/$S.keys.tsv.gz \
    final_all_vcf $U/GATKALL_$P/$S.vcf.gz final_pass_keys $CONC_OUT/prep/GATK_$P/$S.keys.tsv.gz > $R/inputs.tsv

# 3. the status of each raw record in the final file, and the subsets PASS + VQSR, + ExcessHet, + other filters
python3 $HERE/lib/reason_subsets.py --raw-vcf $U/GATKRAW_$P/$S.vcf.gz --raw-keys $U/GATKRAW_$P/$S.keys.tsv.gz \
    --final-all-vcf $U/GATKALL_$P/$S.vcf.gz --final-pass-keys $CONC_OUT/prep/GATK_$P/$S.keys.tsv.gz --out-dir $R/subsets \
    --require-pass-equal

# 4. each subset (baseline) against the PPR calls of the individual, and the classes of the records
for L in $(awk -F'\t' '$6 == "yes" {print $1}' $R/subsets/levels.tsv); do
    rtg vcfeval -b $R/subsets/$L.vcf.gz -c $CONC_OUT/prep/PPR_$P/$S.vcf.gz -t $SDF -o $R/vcfeval_$L \
        --output-mode=ga4gh --ref-overlap -T $THREADS
    touch $R/vcfeval_$L.DONE        # 04_classify_count.py reads only finished runs
    python3 $CONC/04_classify_count.py --a-keys $R/subsets/$L.keys.tsv.gz --b-keys $CONC_OUT/prep/PPR_$P/$S.keys.tsv.gz \
        --vcfeval $R/vcfeval_$L --out-prefix $R/classify_$L
done
