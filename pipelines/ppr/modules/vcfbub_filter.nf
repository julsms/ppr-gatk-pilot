process VCFBUB_FILTER {
    tag "$sample_name"

    // vcfbub 0.1.2 from biocontainers.
    container 'quay.io/biocontainers/vcfbub:0.1.2--hc1c3326_1'

    input:
    tuple val(sample_name), path(vcf)

    output:
    tuple val(sample_name), path("${sample_name}_filtered.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    # The vcfbub biocontainer has no bgzip/tabix -- and vcfbub reads plain VCF
    # fine (verified), so feed it the vg-deconstruct output directly.
    vcfbub -l 0 -a 100000 -i ${vcf} > temp_filtered.vcf

    # Fix chromosome names: remove GRCh38#0# prefix
    sed 's/GRCh38#0#//g' temp_filtered.vcf > temp_fixed_chr.vcf

    # Step 1: remove CHM13 from the graph VCF (column 10 = CHM13).
    # BusyBox cut (vcfbub biocontainer) has no --complement; this awk is
    # byte-identical to GNU 'cut --complement -f10' (## lines pass through).
    awk 'BEGIN{FS=OFS="\\t"} NF<10{print;next} {s=\$1; for(i=2;i<=NF;i++) if(i!=10) s=s OFS \$i; print s}' temp_fixed_chr.vcf > temp_no_chm13.vcf

    # Fix chromosome names to match FASTA format
    sed -e 's/chrM/chrMT/g' \\
        -e 's/chr\\([0-9XY]\\+\\)_GL\\([0-9]\\+\\)v\\([0-9]\\+\\)_random/chrGL\\2.\\3/g' \\
        -e 's/chr\\([0-9XY]\\+\\)_KI\\([0-9]\\+\\)v\\([0-9]\\+\\)_random/chrKI\\2.\\3/g' \\
        -e 's/chrUn_GL\\([0-9]\\+\\)v\\([0-9]\\+\\)/chrGL\\1.\\2/g' \\
        -e 's/chrUn_KI\\([0-9]\\+\\)v\\([0-9]\\+\\)/chrKI\\1.\\2/g' \\
        temp_no_chm13.vcf > temp_fixed_names.vcf

    # Filter to keep only chromosomes that exist in FASTA
    awk '
    BEGIN {
        valid_chrs["chr1"] = 1; valid_chrs["chr2"] = 1; valid_chrs["chr3"] = 1;
        valid_chrs["chr4"] = 1; valid_chrs["chr5"] = 1; valid_chrs["chr6"] = 1;
        valid_chrs["chr7"] = 1; valid_chrs["chr8"] = 1; valid_chrs["chr9"] = 1;
        valid_chrs["chr10"] = 1; valid_chrs["chr11"] = 1; valid_chrs["chr12"] = 1;
        valid_chrs["chr13"] = 1; valid_chrs["chr14"] = 1; valid_chrs["chr15"] = 1;
        valid_chrs["chr16"] = 1; valid_chrs["chr17"] = 1; valid_chrs["chr18"] = 1;
        valid_chrs["chr19"] = 1; valid_chrs["chr20"] = 1; valid_chrs["chr21"] = 1;
        valid_chrs["chr22"] = 1; valid_chrs["chrX"] = 1; valid_chrs["chrY"] = 1;
        valid_chrs["chrMT"] = 1;

        valid_chrs["chrGL000008.2"] = 1; valid_chrs["chrGL000009.2"] = 1;
        valid_chrs["chrGL000194.1"] = 1; valid_chrs["chrGL000205.2"] = 1;
        valid_chrs["chrGL000214.1"] = 1; valid_chrs["chrGL000216.2"] = 1;
        valid_chrs["chrGL000218.1"] = 1; valid_chrs["chrGL000220.1"] = 1;
        valid_chrs["chrGL000224.1"] = 1; valid_chrs["chrGL000225.1"] = 1;
        valid_chrs["chrKI270435.1"] = 1; valid_chrs["chrKI270438.1"] = 1;
        valid_chrs["chrKI270442.1"] = 1; valid_chrs["chrKI270741.1"] = 1;
        valid_chrs["chrKI270742.1"] = 1; valid_chrs["chrKI270743.1"] = 1;
        valid_chrs["chrKI270744.1"] = 1; valid_chrs["chrKI270746.1"] = 1;
        valid_chrs["chrKI270748.1"] = 1; valid_chrs["chrKI270749.1"] = 1;
        valid_chrs["chrKI270751.1"] = 1;
    }
    /^#/ {print; next}
    valid_chrs[\$1] {print}
    ' temp_fixed_names.vcf > ${sample_name}_filtered.vcf

    # Clean up transient intermediates to save work-dir space

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfbub: \$(vcfbub --version 2>&1 | head -1)
    END_VERSIONS
    """
}
