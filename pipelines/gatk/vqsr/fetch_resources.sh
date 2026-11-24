#!/bin/bash
# Downloads the Broad hg38 VQSR resources of WARP JointGenotyping v1.7.3 (WGS) and renames their contigs to the
# Ensembl names of the GATK reference (Homo_sapiens.GRCh38.dna.primary_assembly).
# Usage: bash fetch_resources.sh <reference .dict> <output folder>
set -euo pipefail

DICT=$(readlink -f "$1")
OUT=$2
HERE=$(cd "$(dirname "$0")" && pwd)
GCS=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
mkdir -p $OUT/original
cd $OUT

# Download, then check the size of every object and the md5 where GCS gives one (dbSNP has none;
# its crc32c in resources.tsv can be checked with gsutil hash -c)
tail -n +2 $HERE/resources.tsv | cut -f2 | while read obj; do
    curl -fsSL -o original/$obj $GCS/$obj
done
tail -n +2 $HERE/resources.tsv | while IFS=$'\t' read id obj size md5 crc renamed; do
    [[ $(stat -c %s original/$obj) == "$size" ]] || { echo "size of $obj differs from resources.tsv"; exit 1; }
done
tail -n +2 $HERE/resources.tsv | awk -F'\t' '$4 != "-" {print $4 "  original/" $2}' | md5sum -c -

# chrN to N, chrM to MT, unplaced scaffolds to their accession; records on contigs that are not in the
# reference (alt, HLA, decoy, EBV) are dropped; records sorted in the order of the .dict
tail -n +2 $HERE/resources.tsv | awk -F'\t' '$6 != "-" {print $1, $2, $6}' | while read id obj renamed; do
    gzip -cdf original/$obj \
        | awk -v DICT=$DICT -v STATS=$id.rename_stats.tsv \
              -v NOTE="\"UCSC to Ensembl contig names of $(basename $DICT), contigs and order of that dict, records on other contigs dropped; source $GCS/$obj\"" \
              -f $HERE/rename_contigs.awk \
        | bcftools sort -m 4G -Oz -o $renamed -
    bcftools index -f -t $renamed
done
