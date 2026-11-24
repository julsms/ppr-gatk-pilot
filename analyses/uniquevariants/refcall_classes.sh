#!/bin/bash
# RefCall classes of three unique sets of one individual: is the unique allele a RefCall record (a candidate that DeepVariant
# rejected) of the other call set? GATK vs PPR and DV vs PPR against the PPR calls, PPR vs DV against the DV calls.
# Usage: bash refcall_classes.sh <sample> <platform> <PPR DeepVariant VCF> <DV DeepVariant VCF>    (after extract_unique.sh)
set -euo pipefail

S=$1
P=$2
PPR_VCF=$3                          # <sample>.deepvariant.vcf.gz of the PPR pipeline, RefCall records included
DV_VCF=$4                           # <sample>.deepvariant.vcf.gz of BWA-MEM + DeepVariant
REF=${REF:-GRCh38_chr.fa}           # GRCh38 with chr names (+ .fai), as in the concordance step
UV_OUT=${UV_OUT:-uniquevariants}
CONC=$(cd "$(dirname "$0")/../concordance" && pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
D=$UV_OUT/refcall
# the name of the individual in one call set (column <platform>_ppr or <platform>_dv of the concordance sample map)
name() { awk -F'\t' -v s=$S -v c=$1 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $h["patient"] == s {print $h[c]}' \
             $CONC/config/sample_map.tsv; }
mkdir -p $D/GATKvsPPR_$P $D/DVvsPPR_$P $D/PPRvsDV_$P

# RefCall records of the PPR calls (chr names), split, trimmed and left-aligned as the call sets
bcftools view -s $(name ${P}_ppr) -f RefCall -Ou $PPR_VCF \
    | bcftools annotate -x INFO,^FORMAT/GT -Ou - \
    | bcftools norm -f $REF --check-ref e -m -both -Ou - \
    | bcftools view -T $CONC/config/analysis_space.bed -e 'ALT="*"' -Ou - \
    | bcftools norm -d exact -Ou - \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' - > $D/GATKvsPPR_$P/$S.refcall_keys.tsv

# RefCall records of the DV calls (Ensembl names, renamed to chr names)
bcftools view -s $(name ${P}_dv) -f RefCall -Ou $DV_VCF \
    | bcftools annotate --rename-chrs $CONC/config/ensembl_to_chr.map -x INFO,^FORMAT/GT -Ou - \
    | bcftools norm -f $REF --check-ref e -m -both -Ou - \
    | bcftools view -T $CONC/config/analysis_space.bed -e 'ALT="*"' -Ou - \
    | bcftools norm -d exact -Ou - \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' - > $D/PPRvsDV_$P/$S.refcall_keys.tsv

# classes: i_refcall (the key is a RefCall record of the other set) or ii_none
for set_keys in GATKvsPPR:GATKvsPPR DVvsPPR:GATKvsPPR PPRvsDV:PPRvsDV; do
    U=${set_keys%:*}_$P
    K=$D/${set_keys#*:}_$P/$S.refcall_keys.tsv
    python3 $HERE/lib/refcall_classes.py --unique-keys $UV_OUT/unique/$U/$S.keys.tsv.gz --refcall-keys $K --uset $U --patient $S \
        --out-prefix $D/$U/$S
    touch $D/$U/$S.DONE             # 05_tables.py reads only the individuals with a .DONE file
done
