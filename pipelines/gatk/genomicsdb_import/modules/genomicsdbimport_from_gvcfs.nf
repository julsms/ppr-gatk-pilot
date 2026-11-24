// One new GenomicsDB store for one interval, from the gVCFs of all samples. The GenomicsDBImport options are those of the
// per-sample runs of the cohort (illumina: -Xmx54g with 60 GC threads; aviti: -Xmx12g, --reader-threads 8, --batch-size 50).
// The store is written as <store>.incomplete, checked, then renamed; an existing store is never updated.
process GENOMICSDBIMPORT_FROM_GVCFS {
    tag "${interval}"

    input:
    tuple val(interval), val(store), val(sample_ids), path(gvcfs), path(tbis)

    output:
    tuple val(interval), path("${store}.check.tsv"), emit: check
    tuple val(interval), path("${store}"), emit: store
    path "versions.yml", emit: versions

    script:
    // absolute path: a relative --genomicsdb_dir is resolved against the launch folder
    def store_path = "${file(params.genomicsdb_dir)}/${store}"
    def safe = interval.replaceAll('[^A-Za-z0-9]', '_')
    // Staged gVCF of each sample: <sample>_<interval with each non-alphanumeric character as _>.g.vcf.gz
    def map_args = sample_ids.collect { s -> "'${s}' '${s}_${safe}.g.vcf.gz'" }.join(' ')
    def import_cmd = ''
    if (params.platform == 'illumina') {
        def v_args = sample_ids.collect { s -> "-V ${s}_${safe}.g.vcf.gz" }.join(' \\\n        ')
        import_cmd = [
            "mkdir -p ${params.tmp_dir}",
            "    gatk --java-options \"-XX:+UseParallelGC -XX:ParallelGCThreads=60 -Xmx54g -Djava.io.tmpdir=${params.tmp_dir}\" \\",
            "        GenomicsDBImport \\",
            "        ${v_args} \\",
            "        --genomicsdb-workspace-path \"${store_path}.incomplete\" \\",
            "        -L \"${interval}\" \\",
            "        --tmp-dir ${params.tmp_dir}"
        ].join('\n')
    }
    else {
        import_cmd = [
            "gatk --java-options \"-Xmx12g\" GenomicsDBImport \\",
            "        --genomicsdb-workspace-path \"${store_path}.incomplete\" \\",
            "        --sample-name-map gvcf_map.txt \\",
            "        -L \"${interval}\" \\",
            "        --reader-threads 8 \\",
            "        --batch-size 50"
        ].join('\n')
    }
    """
    store="${store_path}"
    if [[ -e "\$store" || -e "\$store.incomplete" ]]; then
        echo "ERROR: \$store or \$store.incomplete exists. This step never updates a store." >&2
        exit 1
    fi

    # Sample name and staged gVCF of each sample (Aviti reads this map; Illumina gives the same files with -V)
    printf '%s\\t%s\\n' ${map_args} > gvcf_map.txt

    # Each name in the sheet must be the one sample of its gVCF header (the header line only)
    bad=0
    while IFS=\$'\\t' read -r name gvcf; do
        found=\$( (gzip -cd "\$gvcf" 2>/dev/null || true) | awk -F'\\t' '/^#CHROM/ { print (NF == 10 ? \$10 : "not_one_sample"); exit }')
        if [[ "\$found" != "\$name" ]]; then
            echo "ERROR: sample sheet name \$name, gVCF header name \${found:-none} (\$gvcf)" >&2
            bad=1
        fi
    done < gvcf_map.txt
    [[ \$bad -eq 0 ]] || exit 1

    ${import_cmd}

    # The store must hold exactly the samples of the sheet
    cut -f1 gvcf_map.txt | sort > expected_names.txt
    grep -o '"sample_name" *: *"[^"]*"' "\$store.incomplete/callset.json" | sed 's/.*: *"//; s/"\$//' | sort > found_names.txt
    if ! cmp -s expected_names.txt found_names.txt; then
        echo "ERROR: the samples in \$store.incomplete/callset.json differ from the sample sheet:" >&2
        diff expected_names.txt found_names.txt >&2 || true
        exit 1
    fi

    mv "\$store.incomplete" "\$store"
    # A link in the work folder: with -resume, Nextflow runs the task again if the store is gone
    ln -s "\$store" ${store}
    printf '%s\\t%s\\t%s\\t%s\\tOK\\n' "${interval}" "${store}" "\$(wc -l < expected_names.txt)" "\$(wc -l < found_names.txt)" > ${store}.check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def store_path = "${file(params.genomicsdb_dir)}/${store}"
    """
    store="${store_path}"
    if [[ -e "\$store" || -e "\$store.incomplete" ]]; then
        echo "ERROR: \$store or \$store.incomplete exists. This step never updates a store." >&2
        exit 1
    fi
    mkdir -p "\$store"
    printf '%s\\n' ${sample_ids.join(' ')} > "\$store/stub_sample_names.txt"
    ln -s "\$store" ${store}
    printf '%s\\t%s\\t%s\\t%s\\tSTUB\\n' "${interval}" "${store}" "${sample_ids.size()}" "${sample_ids.size()}" > ${store}.check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: stub
    END_VERSIONS
    """
}
