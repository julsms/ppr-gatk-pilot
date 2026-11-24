#!/bin/bash
# k* of the HG002 calls on one platform: copy number of the k-mers of every record of the three full sets in the reads and in
# GRCh38 (meryl-lookup counts, turned into the values of merfin -dump with the haploid peak of GenomeScope 2.0), the flags, and
# the Merfin x GIAB tables split by them. Run after merfin_sample.sh on HG002 (it uses that read database) and tables.sh.
# Usage: bash kstar_hg002.sh <platform> <output folder> <reference fasta> <reference k-mers> <unique-variant output folder>
set -euo pipefail

P=$1
OUT=$2
REF=$3
SEQMERS=$4
UV=$5                                       # hg002/truth/<set>_<platform>.status.tsv.gz of analyses/uniquevariants
GENOMESCOPE=${GENOMESCOPE:-genomescope.R}   # genomescope.R of GenomeScope 2.0 (tag v2.1.0)
THREADS=8
MEM=16                                      # GB
HERE=$(cd "$(dirname "$0")" && pwd)
DB=$OUT/db/${P}_HG002/reads.gt1.meryl
KS=$OUT/hg002/kstar/$P
W=$KS/work
export OMP_NUM_THREADS=$THREADS
mkdir -p $W $KS/logs

# k-mer histogram of the reads, GenomeScope 2.0 fit, haploid peak (kcov when within 0.75-1.33 x the histogram estimate)
meryl histogram $DB 2>&1 | grep -E '^[0-9]+[[:space:]]+[0-9]+$' > $KS/hist.tsv
Rscript $GENOMESCOPE -i $KS/hist.tsv -o $KS/genomescope -k 31 -p 2 --fitted_hist > $KS/logs/genomescope.log 2>&1
python3 $HERE/lib/kstar.py peak --hist $KS/hist.tsv --genomescope $KS/genomescope --out $KS/peak.tsv
PEAK=$(awk -F'\t' '$1 == "haploid_peak" {print $2}' $KS/peak.tsv)

# REF and ALT windows of the records (30 bases on each side), their k-mers, and both databases reduced to these k-mers
keys=$(awk -F'\t' -v p=$P '$1 == p && $2 == "HG002" && $4 == "all" {print $9}' $OUT/sheet/vcfs.tsv | paste -sd,)
python3 $HERE/lib/kstar.py windows --keys $keys --ref $REF --k 31 --fasta $W/windows.fa --index $W/windows.tsv.gz \
    > $KS/logs/windows.log 2>&1
meryl count k=31 threads=$THREADS memory=$MEM $W/windows.fa output $W/windows.meryl
meryl intersect $DB $W/windows.meryl output $W/reads_sub.meryl
meryl intersect $SEQMERS $W/windows.meryl output $W/ref_sub.meryl

# Count of every window k-mer in the reads and in GRCh38, both at once
meryl-lookup -wig-count -sequence $W/windows.fa -mers $W/reads_sub.meryl -output $W/lookup_reads.wig.gz \
    -threads 4 -memory $MEM > $KS/logs/lookup_reads.log 2>&1 &
p1=$!
meryl-lookup -wig-count -sequence $W/windows.fa -mers $W/ref_sub.meryl -output $W/lookup_ref.wig.gz \
    -threads 4 -memory $MEM > $KS/logs/lookup_ref.log 2>&1 &
p2=$!
wait $p1
wait $p2

# merfin -dump of every n-th record (about 2,000): kstar.py kmers checks that the counts give the same lines
n=$(( $(zcat $W/windows.tsv.gz | wc -l) - 1 ))
every=$(( n > 2000 ? n / 2000 : 1 ))
awk -v n=$every '/^>/ {split(substr($1, 2), a, "_"); keep = (a[1] % n == 0)} keep' $W/windows.fa > $W/spot.fa
merfin -dump -sequence $W/spot.fa -seqmers $W/ref_sub.meryl -readmers $W/reads_sub.meryl -peak $PEAK \
    -output $W/spot_dump.tsv.gz -threads $THREADS -memory $MEM 2>&1 > /dev/null \
    | grep -v -E '^(Processing|Output) sequence |^[0-9]+_[AR][[:space:]]' > $KS/logs/merfin_spot.log

# Copy numbers per record, then the flags and the tables (HG002 is male: one copy on chrY and on chrX outside XPAR)
python3 $HERE/lib/kstar.py kmers --reads-wig $W/lookup_reads.wig.gz --ref-wig $W/lookup_ref.wig.gz --peak $PEAK \
    --index $W/windows.tsv.gz --spot-fasta $W/spot.fa --spot-dump $W/spot_dump.tsv.gz --out $KS/kmer_stats.tsv.gz \
    > $KS/logs/kmers.log 2>&1
python3 $HERE/lib/kstar.py tables --sheet-dir $OUT/sheet --results $OUT/results --state $OUT/sample_state.tsv \
    --sets $HERE/config/merfin_sets.tsv --uv-out $UV --platform $P --sample HG002 --kstats $KS/kmer_stats.tsv.gz \
    --out $KS --kstar-min 1 --min-kmers 5 --haploid-xy

# merfin_table.py reads only the k* folders with this file
printf 'peak\t%s\n' $PEAK > $KS/KSTAR_DONE

# Copy numbers of the supported GIAB true positives in single-copy context (the medians of the paper) and the 2x2 against tables.sh
python3 $HERE/lib/kstar.py verify --kstar $KS --hg002 $OUT/hg002 --platform $P --tp-max 0.15 > $KS/verify.txt 2>&1
