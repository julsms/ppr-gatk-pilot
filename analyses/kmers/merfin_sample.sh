#!/bin/bash
# K-mer support of the call sets of one sample (one person on one platform): read k-mer database (meryl, k = 31, k-mers seen
# more than once), merfin -filter on each call set of the sheet, the status of every record and the counts by region and type.
# Usage: bash merfin_sample.sh <platform> <sample> <output folder> <reference fasta> <reference k-mers>
set -euo pipefail

P=$1
S=$2
OUT=$3                      # holds sheet/samples.tsv and sheet/vcfs.tsv
REF=$4                      # GRCh38 with chr names, with its .fai
SEQMERS=$5                  # output of ref_kmers.sh
THREADS=${THREADS:-16}      # meryl count and merfin
MEM=${MEM:-100}             # GB for meryl count
HERE=$(cd "$(dirname "$0")" && pwd)
RES=$OUT/results/$P/$S
DB=$OUT/db/${P}_$S
read -r R1 R2 < <(awk -F'\t' -v p=$P -v s=$S '$1 == p && $2 == s {print $4, $5}' $OUT/sheet/samples.tsv)
mkdir -p $RES $DB

# Read k-mer database: one count of R1 and R2, then the k-mers seen more than once
meryl count k=31 threads=$THREADS memory=$MEM $R1 $R2 output $DB/reads.meryl
meryl greater-than 1 $DB/reads.meryl output $DB/reads.gt1.meryl
rm -r $DB/reads.meryl

# merfin -filter on each call set (defaults: no -peak, -prob or -comb), then the status of every record (1 = kept by merfin)
awk -F'\t' -v p=$P -v s=$S '$1 == p && $2 == s && $7 == "yes" {print $3, $8}' $OUT/sheet/vcfs.tsv | while read -r id vcf; do
    merfin -filter -sequence $REF -seqmers $SEQMERS -readmers $DB/reads.gt1.meryl -vcf $vcf -output $RES/$id \
        -threads $THREADS > $RES/$id.merfin.log 2>&1
    python3 $HERE/lib/merfin_status.py set --input $vcf --validated $RES/$id.filter.vcf --out $RES/$id.status.tsv.gz \
        --bcftools bcftools
    bgzip -f -c $RES/$id.filter.vcf > $RES/$id.validated.vcf.gz
    tabix -f -p vcf $RES/$id.validated.vcf.gz
    rm $RES/$id.filter.vcf
done

# Counts by region and type: direct (each set alone) and context (the unique records in the run of their full set)
python3 $HERE/lib/merfin_status.py sample --vcfs $OUT/sheet/vcfs.tsv --platform $P --sample $S --results $RES \
    --out $RES/sample_counts.tsv

# The table scripts and lib/kstar.py read only the samples recorded here (SAMPLE_DONE and sample_state.tsv)
fp=$(md5sum < $RES/sample_counts.tsv | cut -c1-32)
printf 'fingerprint\t%s\n' $fp > $RES/SAMPLE_DONE
[[ -s $OUT/sample_state.tsv ]] || printf 'platform\tsample\texit\tstage\tdate\tsample_done_fp\n' > $OUT/sample_state.tsv
printf '%s\t%s\t0\tall\t%s\t%s\n' $P $S "$(date -u '+%F %T UTC')" $fp >> $OUT/sample_state.tsv
