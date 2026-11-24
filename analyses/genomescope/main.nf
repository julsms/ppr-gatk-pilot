#!/usr/bin/env nextflow

/*
 * GenomeScope 2.0 on the reads of one sample: rasusa subsampling (optional), FastK (k = 21), Histex -G, GenomeScope 2.0.
 * Usage: nextflow run main.nf --sample_id X --fastq_r1 R1 --fastq_r2 R2 [--subsample true --target_bases N] --output_dir DIR
 */

// ============================================================================
// PARAMETERS
// ============================================================================

params.sample_id    = null
params.fastq_r1     = null
params.fastq_r2     = null
params.output_dir   = "./genomescope_results"
params.k            = 21
params.ploidy       = 2
params.threads      = 35
params.memory_gb    = 70

// Subsampling (set by run_genomescope.sh from the depth of both platforms)
params.subsample    = false
params.target_bases = null  // total bases target for this sample

// ============================================================================
// PROCESSES
// ============================================================================


process SUBSAMPLE_READS {
    tag "${sample_id}"
    cpus 60
    memory '90 GB'
    time '12h'
    container 'community.wave.seqera.io/library/rasusa_fastk_seqtk:c8764ae5239ee72b'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_sub_R1.fq.gz"), path("${sample_id}_sub_R2.fq.gz")

    script:
    """
    echo "=== Subsampling ${sample_id} to ${params.target_bases} total bases with rasusa ==="

    rasusa reads \
        --bases ${params.target_bases} \
        --seed 42 \
        -o ${sample_id}_sub_R1.fq.gz \
        -o ${sample_id}_sub_R2.fq.gz \
        ${r1} ${r2}

    echo "Subsampling done"
    """
}


process FASTK_COUNT {
    tag "${sample_id}"
    cpus params.threads
    memory "${params.memory_gb} GB"
    time '8h'
    container 'community.wave.seqera.io/library/rasusa_fastk_seqtk:c8764ae5239ee72b'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.hist")

    script:
    """
    export TMPDIR=\$(pwd)

    echo "=== FastK (k=${params.k}, threads=${task.cpus}) for ${sample_id} ==="
    echo "Working directory: \$(pwd)"
    df -h .

    FastK -v -T${task.cpus} -k${params.k} -M${params.memory_gb} -P. -N${sample_id} ${r1} ${r2}

    # Verify histogram was created
    if [[ ! -f "${sample_id}.hist" ]]; then
        echo "ERROR: FastK did not produce ${sample_id}.hist"
        ls -lha
        exit 1
    fi

    echo "FastK done. Histogram:"
    ls -lha ${sample_id}.hist
    """
}

process HISTEX_EXPORT {
    tag "${sample_id}"
    cpus 2
    memory '4 GB'
    time '30m'
    container 'community.wave.seqera.io/library/rasusa_fastk_seqtk:c8764ae5239ee72b'
    publishDir "${params.output_dir}/${sample_id}", mode: 'copy', pattern: "*.histo"

    input:
    tuple val(sample_id), path(hist_file)

    output:
    tuple val(sample_id), path("${sample_id}_k${params.k}.histo")

    script:
    """
    echo "=== Histex -G for ${sample_id} ==="

    # -G: GenomeScope format; the last bin holds all k-mers at or above it.
    Histex -G ${sample_id} > ${sample_id}_k${params.k}.histo

    echo "First 10 lines:"
    head -10 ${sample_id}_k${params.k}.histo
    echo "Last 5 lines:"
    tail -5 ${sample_id}_k${params.k}.histo
    echo "Total lines: \$(wc -l < ${sample_id}_k${params.k}.histo)"
    echo "Histogram exported"
    """
}

process GENOMESCOPE {
    tag "${sample_id}"
    cpus 2
    memory '4 GB'
    time '30m'
    container 'community.wave.seqera.io/library/genomescope2:2.1.0--c39bc08092274c2b'
    publishDir "${params.output_dir}/${sample_id}", mode: 'copy', pattern: "genomescope2_output/**"

    input:
    tuple val(sample_id), path(histogram)

    output:
    tuple val(sample_id), path("genomescope2_output/**")

    script:
    """
    echo "=== GenomeScope 2.0 for ${sample_id} ==="
    genomescope2 \
        -i ${histogram} \
        -o genomescope2_output \
        -k ${params.k} \
        -p ${params.ploidy} \
        -n ${sample_id}

    if [[ -f "genomescope2_output/summary.txt" ]]; then
        echo "=== Summary ==="
        cat genomescope2_output/summary.txt
    fi
    echo "GenomeScope done"
    """
}

// ============================================================================
// WORKFLOW
// ============================================================================

workflow {
    if (!params.sample_id || !params.fastq_r1 || !params.fastq_r2) {
        error "Required: --sample_id, --fastq_r1, --fastq_r2"
    }

    def input_ch = channel.of(
        [params.sample_id, file(params.fastq_r1), file(params.fastq_r2)]
    )

    def reads_ch
    if (params.subsample && params.target_bases) {
        reads_ch = SUBSAMPLE_READS(input_ch)
    } else {
        reads_ch = input_ch
    }

    def fastk_out = FASTK_COUNT(reads_ch)
    def histo_ch  = HISTEX_EXPORT(fastk_out)
    GENOMESCOPE(histo_ch)
}
