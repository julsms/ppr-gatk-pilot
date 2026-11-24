process GENOMICSDBIMPORT {
    tag "${sample_id}:${interval}"

    input:
    tuple val(sample_id), val(interval), path(gvcf), path(tbi)

    output:
    tuple val(sample_id), val(interval), path("genomicsdb_${interval.replaceAll('[^A-Za-z0-9]', '_')}"), emit: genomicsdb
    path "versions.yml", emit: versions

    script:
    def interval_safe = interval.replaceAll('[^A-Za-z0-9]', '_')
    // One store per interval, shared by every sample: a sample is added to the store if it exists
    def persistent_db = "${params.genomicsdb_dir}/genomicsdb_${interval_safe}"
    """
    echo "=== GenomicsDBImport ==="
    echo "Sample: ${sample_id}"
    echo "Interval: ${interval}"
    echo "GVCF: ${gvcf}"
    echo "Persistent DB: ${persistent_db}"
    echo "CPUs: ${task.cpus}"
    echo "Memory: ${task.memory.toGiga()}GB"
    echo "========================"
    
    # Create a map file for the GVCF
    echo -e "${sample_id}\t${gvcf}" > gvcf_map.txt
    
    echo "GVCF map file contents:"
    cat gvcf_map.txt
    
    # Check if GenomicsDB workspace already exists
    if [ -d "${persistent_db}" ]; then
        echo "GenomicsDB workspace already exists at ${persistent_db}"
        echo "UPDATING existing workspace with new sample: ${sample_id}"
        
        # Update existing GenomicsDB workspace
        gatk --java-options "-Xmx${task.memory.toGiga() - 4}g" GenomicsDBImport \\
            --genomicsdb-update-workspace-path ${persistent_db} \\
            --sample-name-map gvcf_map.txt \\
            -L "${interval}" \\
            --reader-threads ${task.cpus} \\
            --batch-size 50
        
        # Copy updated workspace to work directory for output
        cp -r ${persistent_db} genomicsdb_${interval_safe}
        
        echo "GenomicsDB workspace UPDATED with sample ${sample_id}"
    else
        echo "Creating NEW GenomicsDB workspace at ${persistent_db}"
        
        # Create new GenomicsDB workspace
        gatk --java-options "-Xmx${task.memory.toGiga() - 4}g" GenomicsDBImport \\
            --genomicsdb-workspace-path genomicsdb_${interval_safe} \\
            --sample-name-map gvcf_map.txt \\
            -L "${interval}" \\
            --reader-threads ${task.cpus} \\
            --batch-size 50
        
        # Copy to persistent storage
        mkdir -p \$(dirname ${persistent_db})
        cp -r genomicsdb_${interval_safe} ${persistent_db}
        
        echo "GenomicsDB workspace CREATED and saved to persistent storage"
    fi
    
    echo "GenomicsDBImport completed successfully"
    echo "GenomicsDB directory: genomicsdb_${interval_safe}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'Version:\\s*\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """
}
