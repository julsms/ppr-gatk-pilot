#!/bin/bash
# The GATK status of the PPR-only alleles of every individual on one platform: alleles that no raw call of the individual explains
# are looked up in the raw joint calls of the cohort (0/0 or missing in this individual, or not called at all), then the classes.
# Usage: bash gatk_status_cohort.sh <platform> <raw joint VCF>    (after gatk_status.sh for every individual)
set -euo pipefail

P=$1
RAW=$2                              # merged_raw_variants.vcf.gz of the GATK joint calling
REF=${REF:-GRCh38_chr.fa}           # GRCh38 with chr names (+ .fai), as in the concordance step
UV_OUT=${UV_OUT:-uniquevariants}
CONC=$(cd "$(dirname "$0")/../concordance" && pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
RD=$UV_OUT/reasons/new_$P
# a column of the concordance sample map, for the individuals that are not excluded
col() { awk -F'\t' -v c=$1 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $h["exclude"] != "yes" {print $h[c]}' \
            $CONC/config/sample_map.tsv; }
SAMPLES=$(col patient)

# the candidate alleles of every individual, and windows of 1 kb on both sides of them
for S in $SAMPLES; do printf '%s\t%s\n' $S $RD/$S; done > $RD/states.manifest.tsv
python3 $HERE/lib/reason_states.py collect --manifest $RD/states.manifest.tsv --out-prefix $RD/states

# their genotypes in the raw joint calls of the cohort, split and left-aligned as the call sets (one pass over the file)
paste <(col ${P}_gatk) <(col patient) > $RD/states.labels.tsv
bcftools view -T $RD/states.targets.bed --targets-overlap 1 -s $(cut -f1 $RD/states.labels.tsv | paste -sd,) -Ou $RAW \
    | bcftools annotate --rename-chrs $CONC/config/ensembl_to_chr.map -x INFO,^FORMAT/GT -Ou - \
    | bcftools norm -f $REF --check-ref e -m -both -Ou - \
    | bcftools view -e 'ALT="*"' -Ou - \
    | bcftools query -H -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n' - \
    | python3 $HERE/lib/reason_states.py states --candidates $RD/states.candidates.tsv.gz --label-map $RD/states.labels.tsv \
        --out $RD/states.tsv.gz

# the class of every PPR-only allele of every individual
for S in $SAMPLES; do
    python3 $HERE/lib/reason_classes.py --eval-dir $RD/$S --states $RD/states.tsv.gz --patient $S --mode new --platform $P \
        --out-prefix $RD/$S
    touch $RD/$S.DONE               # 05_tables.py reads only the individuals with a .DONE file
done
