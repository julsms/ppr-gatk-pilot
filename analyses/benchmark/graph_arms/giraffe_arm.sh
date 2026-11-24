#!/bin/bash
# One giraffe set-up of the graph decomposition for one sample: vg giraffe on the given graph, vg stats -a on the alignments,
# vg surject onto the 24 GRCh38 paths with the read group and sort of the PPR pipeline, then DeepVariant 1.9.0 (PPR settings).
# Usage: bash giraffe_arm.sh <arm folder> <sample> <R1.fq.gz> <R2.fq.gz> <graph .gbz> <.dist> <.min> <GRCh38 FASTA, chr names>
set -euo pipefail

A=$1                    # for example Illumina_giraffe_HPRC_d9_Deepvariant, made in the current folder
S=$2                    # for example HG002Illumina_SRR14724532
R1=$(readlink -f $3)
R2=$(readlink -f $4)
GBZ=$(readlink -f $5)
DIST=$(readlink -f $6)
MIN=$(readlink -f $7)
REF=$(readlink -f $8)   # with its .fai
THREADS=${THREADS:-64}
VGSTATS_THREADS=${VGSTATS_THREADS:-4}
SORT_THREADS=16
SORT_MEM=4G
DV_SHARDS=40
VG=quay.io/vgteam/vg:v1.52.0@sha256:795887f8381623d0242743b5fbb813e0fcc1b0fdd8dae04f335e0798f1e454e5
VGS=ppr/vg-samtools:1.52.0-1.17   # docker build -t ppr/vg-samtools:1.52.0-1.17 <release>/pipelines/ppr/docker/vg-samtools
DV=google/deepvariant:1.9.0-gpu@sha256:e0c8734b8700d945e3ee78d609acb90548f829c874596ffca436af8cf379f87a
P=""
for c in $(seq 1 22) X Y; do P+=" -p GRCh38#0#chr$c"; done
mkdir -p $A/work/sorttmp $A/mapping $A/out
# each step runs in its own container; the current folder and the folders of the inputs keep their paths there
MOUNTS=$(for d in $PWD $(dirname $R1 $R2 $GBZ $DIST $MIN $REF); do echo "-v $d:$d"; done | sort -u)
dock() { docker run --rm -u $(id -u):$(id -g) -e HOME=/tmp $MOUNTS -w $PWD "$@"; }

# giraffe, default preset, the distance and minimizer indexes of the graph given explicitly
dock $VG bash -c "vg giraffe -p -t $THREADS -Z $GBZ -d $DIST -m $MIN -f $R1 -f $R2 > $A/work/$S.gam 2> $A/work/${S}_align.err"

# graph-space alignment metrics (Figure 2e, Supplementary Figure S7)
dock $VG bash -c "vg stats -p $VGSTATS_THREADS -a $A/work/$S.gam > $A/mapping/vg_stats_a.txt"

# surject from the GBZ onto GRCh38, contigs renamed to chr1..chrY, read group, sort: one pipe in one container, which stops at
# the first failed command (errexit, nounset and pipefail, as in the run)
dock $VGS bash -c "set -euo pipefail
vg surject -x $GBZ -b $A/work/$S.gam -t $THREADS --prune-low-cplx -i --max-frag-len 3000 $P 2> $A/out/${S}_surject.log \
    | samtools reheader -c \"sed 's/SN:GRCh38#0#chr/SN:chr/'\" - \
    | samtools addreplacerg -u -r '@RG\tID:$S\tSM:$S\tPL:ILLUMINA\tLB:lib1\tPU:unit1' -o - - \
    | samtools sort -@ $SORT_THREADS -m $SORT_MEM -T $A/work/sorttmp/$S.sorttmp -o $A/out/${S}_chr_named_rg.bam -
samtools index -@ $THREADS $A/out/${S}_chr_named_rg.bam
samtools flagstat -@ $THREADS $A/out/${S}_chr_named_rg.bam > $A/out/$S.flagstat.txt"

# DeepVariant 1.9.0 on the GPU, WGS model, the make_examples settings of the PPR pipeline
dock --gpus all $DV /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=$REF --reads=$A/out/${S}_chr_named_rg.bam \
    --output_vcf=$A/out/$S.deepvariant.vcf.gz --output_gvcf=$A/out/$S.deepvariant.g.vcf.gz --num_shards=$DV_SHARDS \
    --intermediate_results_dir=$A/work/dv_tmp \
    --make_examples_extra_args="min_mapping_quality=1,keep_legacy_allele_counter_behavior=true,normalize_reads=true"
rm -rf $A/work/dv_tmp
