# Joint genotyping, Illumina

GenotypeGVCFs on each of the 66 GenomicsDB stores of `../genomicsdb_import` (`--platform illumina`), then Picard
MergeVcfs. The output `results_joint_analysis/merged_vcf/merged_raw_variants.vcf.gz` is the raw joint VCF of the paper:
the input of `../vqsr`.

GenotypeGVCFs runs with `-G StandardAnnotation -G AS_StandardAnnotation --use-new-qual-calculator
--annotate-with-num-discovered-alleles` and `-A` Coverage, QualByDepth, FisherStrand, StrandOddsRatio, MappingQuality,
MappingQualityRankSumTest, ReadPosRankSumTest, RMSMappingQuality and InbreedingCoeff, with no `-L`.

The VQSR of the paper is `../vqsr`. The runs also had, after MergeVcfs, an older VQSR on split SNP and INDEL files and
an annotation step; their outputs were not used, and these steps are left out here.

## Run

    nextflow run pipeline2.nf -params-file params.yaml \
        --genomicsdb_dir <folder with the 66 stores> --samplename_joint <name>

`params.yaml` gives the paths named in `nextflow.config`: the reference with `.fai` and `.dict`, `picard` and `tmp_dir`.
With `-profile docker --data_dir <folder>` the GATK steps run in the GATK image, and all paths must be under that folder.
The outputs in `results_joint_analysis/` are links into `work/`: copy them with `cp -L`.

## As run

The joint genotyping of HG002 and of HG001 with the 24 cohort genomes ran this code: the GenotypeGVCFs and MergeVcfs task
scripts of each execution report match it. HG001 alone (one sample) ran the same code. In the runs the scratch folder was a
fixed path of the server; here it is the param `tmp_dir`.

Tools: Nextflow 25.04.7, GATK 4.5.0.0, Picard 2.26.0.
