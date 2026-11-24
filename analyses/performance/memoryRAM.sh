#!/bin/bash
# Peak memory per sample from the Nextflow traces: the largest peak_rss of one task and the sum over tasks, without the QC steps
# of GATK and the structural-variant steps of PPR. Output: memory_usage_results.csv.
# Usage: bash memoryRAM.sh [folder with PPR/<platform>/ and GATK/<platform>/individual/ traces]

BASE_DIR=${1:-.}
cd $BASE_DIR

echo "=== CALCULATE PEAK MEMORY USAGE PER SAMPLE ==="

# Function that calculates the peak memory from a trace file
calculate_peak_memory() {
    local trace_file="$1"
    local pipeline_type="$2"
    
    # Find the columns for the pipeline type
    if [[ "$pipeline_type" == "gatk_aviti" ]]; then
        process_col=4
        peak_rss_col=26  # Column 26 for GATK Aviti (peak_rss; column 24 is rss)
    elif [[ "$pipeline_type" == "gatk_illumina" ]]; then
        process_col=4
        peak_rss_col=11  # Column 11 for GATK Illumina (peak_rss)
    elif [[ "$pipeline_type" == "ppr" ]]; then
        process_col=4
        peak_rss_col=11  # Column 11 for PPR (peak_rss)
    fi
    
    awk -F'\t' -v proc_col="$process_col" -v rss_col="$peak_rss_col" '
    BEGIN { max_peak_rss = 0; total_peak_rss = 0 }
    
    # Skip the header and the QC/SV processes
    NR > 1 {
        process_name = $proc_col
        
        # Get only the process name
        if (match(process_name, /^([A-Z_]+)/)) {
            process = substr(process_name, RSTART, RLENGTH)
        } else {
            process = process_name
        }
        
        # Skip the QC and SV calling processes
        if (process !~ /^(MOSDEPTH|SAMTOOLS_STATS|PICARD_ALIGNMENT_METRICS|VG_CALL|VG_PACK|VG_DECONSTRUCT|VCFBUB_FILTER|VCFWAVE|PANGENIE_PREP|PANGENIE_INDEX|PANGENIE_GENOTYPE|PANGENIE_BIALLELIC|COMPRESS_INDEX_VCF_PANGENIE)$/) {
            peak_rss = $rss_col
            
            if (peak_rss != "-" && peak_rss != "" && peak_rss != "peak_rss") {
                # Convert to GB
                if (match(peak_rss, /([0-9]+\.?[0-9]*) GB/)) {
                    rss_gb = substr(peak_rss, RSTART, RLENGTH-3) + 0
                    total_peak_rss += rss_gb
                    if (rss_gb > max_peak_rss) {
                        max_peak_rss = rss_gb
                    }
                } else if (match(peak_rss, /([0-9]+\.?[0-9]*) MB/)) {
                    rss_mb = substr(peak_rss, RSTART, RLENGTH-3) + 0
                    rss_gb = rss_mb / 1024
                    total_peak_rss += rss_gb
                    if (rss_gb > max_peak_rss) {
                        max_peak_rss = rss_gb
                    }
                } else if (match(peak_rss, /([0-9]+\.?[0-9]*)GB/)) {
                    # For the case without a space
                    rss_gb = substr(peak_rss, RSTART, RLENGTH-2) + 0
                    total_peak_rss += rss_gb
                    if (rss_gb > max_peak_rss) {
                        max_peak_rss = rss_gb
                    }
                } else if (match(peak_rss, /([0-9]+\.?[0-9]*)MB/)) {
                    # For the case without a space
                    rss_mb = substr(peak_rss, RSTART, RLENGTH-2) + 0
                    rss_gb = rss_mb / 1024
                    total_peak_rss += rss_gb
                    if (rss_gb > max_peak_rss) {
                        max_peak_rss = rss_gb
                    }
                }
            }
        }
    }
    
    END { 
        print max_peak_rss "," total_peak_rss 
    }
    ' "$trace_file"
}

# Create the results file for the memory usage
results_file="memory_usage_results.csv"
echo "Pipeline,Sample,Max_Peak_RSS_GB,Total_Peak_RSS_GB" > $results_file

# Arrays that store the results
declare -a ppr_illumina_max
declare -a ppr_illumina_total
declare -a ppr_aviti_max
declare -a ppr_aviti_total
declare -a gatk_aviti_max
declare -a gatk_aviti_total
declare -a gatk_illumina_max
declare -a gatk_illumina_total

echo ""
echo "Processing PPR Illumina samples..."
for trace_file in PPR/illumina/*_trace.txt; do
    if [[ -f "$trace_file" ]]; then
        sample_name=$(basename "$trace_file" _trace.txt)
        memory_result=$(calculate_peak_memory "$trace_file" "ppr")
        
        max_rss=$(echo "$memory_result" | cut -d',' -f1)
        total_rss=$(echo "$memory_result" | cut -d',' -f2)
        
        ppr_illumina_max+=($max_rss)
        ppr_illumina_total+=($total_rss)
        
        echo "  $sample_name: Max=${max_rss}GB, Total=${total_rss}GB"
        echo "PPR_Illumina,$sample_name,$max_rss,$total_rss" >> $results_file
    fi
done

echo ""
echo "Processing PPR Aviti samples..."
for trace_file in PPR/aviti/*_trace.txt; do
    if [[ -f "$trace_file" ]]; then
        sample_name=$(basename "$trace_file" _trace.txt)
        memory_result=$(calculate_peak_memory "$trace_file" "ppr")
        
        max_rss=$(echo "$memory_result" | cut -d',' -f1)
        total_rss=$(echo "$memory_result" | cut -d',' -f2)
        
        ppr_aviti_max+=($max_rss)
        ppr_aviti_total+=($total_rss)
        
        echo "  $sample_name: Max=${max_rss}GB, Total=${total_rss}GB"
        echo "PPR_Aviti,$sample_name,$max_rss,$total_rss" >> $results_file
    fi
done

echo ""
echo "Processing GATK Aviti samples..."
for trace_file in GATK/aviti/individual/*_trace.txt; do
    if [[ -f "$trace_file" ]]; then
        sample_name=$(basename "$trace_file" _trace.txt)
        memory_result=$(calculate_peak_memory "$trace_file" "gatk_aviti")
        
        max_rss=$(echo "$memory_result" | cut -d',' -f1)
        total_rss=$(echo "$memory_result" | cut -d',' -f2)
        
        gatk_aviti_max+=($max_rss)
        gatk_aviti_total+=($total_rss)
        
        echo "  $sample_name: Max=${max_rss}GB, Total=${total_rss}GB"
        echo "GATK_Aviti,$sample_name,$max_rss,$total_rss" >> $results_file
    fi
done

echo ""
echo "Processing GATK Illumina samples..."
for trace_file in GATK/illumina/individual/*_trace.txt; do
    if [[ -f "$trace_file" ]]; then
        sample_name=$(basename "$trace_file" _trace.txt)
        # Skip the problem sample
        if [[ "$sample_name" == "23" ]]; then
            continue
        fi

        memory_result=$(calculate_peak_memory "$trace_file" "gatk_illumina")
        
        max_rss=$(echo "$memory_result" | cut -d',' -f1)
        total_rss=$(echo "$memory_result" | cut -d',' -f2)
        
        gatk_illumina_max+=($max_rss)
        gatk_illumina_total+=($total_rss)
        
        echo "  $sample_name: Max=${max_rss}GB, Total=${total_rss}GB"
        echo "GATK_Illumina,$sample_name,$max_rss,$total_rss" >> $results_file
    fi
done

# Calculate the means
echo ""
echo "=== CALCULATE MEAN MEMORY USAGE ==="

# PPR Illumina
if [[ ${#ppr_illumina_max[@]} -gt 0 ]]; then
    total_max=0
    total_sum=0
    for max in "${ppr_illumina_max[@]}"; do
        total_max=$(echo "$total_max + $max" | bc)
    done
    for sum in "${ppr_illumina_total[@]}"; do
        total_sum=$(echo "$total_sum + $sum" | bc)
    done
    avg_max_ppr_illumina=$(echo "scale=1; $total_max / ${#ppr_illumina_max[@]}" | bc)
    avg_total_ppr_illumina=$(echo "scale=1; $total_sum / ${#ppr_illumina_total[@]}" | bc)
    
    echo "PPR Illumina - Mean max: ${avg_max_ppr_illumina}GB, Mean total: ${avg_total_ppr_illumina}GB"
    echo "PPR_Illumina,AVERAGE,$avg_max_ppr_illumina,$avg_total_ppr_illumina" >> $results_file
fi

# PPR Aviti
if [[ ${#ppr_aviti_max[@]} -gt 0 ]]; then
    total_max=0
    total_sum=0
    for max in "${ppr_aviti_max[@]}"; do
        total_max=$(echo "$total_max + $max" | bc)
    done
    for sum in "${ppr_aviti_total[@]}"; do
        total_sum=$(echo "$total_sum + $sum" | bc)
    done
    avg_max_ppr_aviti=$(echo "scale=1; $total_max / ${#ppr_aviti_max[@]}" | bc)
    avg_total_ppr_aviti=$(echo "scale=1; $total_sum / ${#ppr_aviti_total[@]}" | bc)
    
    echo "PPR Aviti - Mean max: ${avg_max_ppr_aviti}GB, Mean total: ${avg_total_ppr_aviti}GB"
    echo "PPR_Aviti,AVERAGE,$avg_max_ppr_aviti,$avg_total_ppr_aviti" >> $results_file
fi

# GATK Aviti
if [[ ${#gatk_aviti_max[@]} -gt 0 ]]; then
    total_max=0
    total_sum=0
    for max in "${gatk_aviti_max[@]}"; do
        total_max=$(echo "$total_max + $max" | bc)
    done
    for sum in "${gatk_aviti_total[@]}"; do
        total_sum=$(echo "$total_sum + $sum" | bc)
    done
    avg_max_gatk_aviti=$(echo "scale=1; $total_max / ${#gatk_aviti_max[@]}" | bc)
    avg_total_gatk_aviti=$(echo "scale=1; $total_sum / ${#gatk_aviti_total[@]}" | bc)
    
    echo "GATK Aviti - Mean max: ${avg_max_gatk_aviti}GB, Mean total: ${avg_total_gatk_aviti}GB"
    echo "GATK_Aviti,AVERAGE,$avg_max_gatk_aviti,$avg_total_gatk_aviti" >> $results_file
fi

# GATK Illumina
if [[ ${#gatk_illumina_max[@]} -gt 0 ]]; then
    total_max=0
    total_sum=0
    for max in "${gatk_illumina_max[@]}"; do
        total_max=$(echo "$total_max + $max" | bc)
    done
    for sum in "${gatk_illumina_total[@]}"; do
        total_sum=$(echo "$total_sum + $sum" | bc)
    done
    avg_max_gatk_illumina=$(echo "scale=1; $total_max / ${#gatk_illumina_max[@]}" | bc)
    avg_total_gatk_illumina=$(echo "scale=1; $total_sum / ${#gatk_illumina_total[@]}" | bc)
    
    echo "GATK Illumina - Mean max: ${avg_max_gatk_illumina}GB, Mean total: ${avg_total_gatk_illumina}GB"
    echo "GATK_Illumina,AVERAGE,$avg_max_gatk_illumina,$avg_total_gatk_illumina" >> $results_file
fi

echo ""
echo "The results were saved to: $results_file"
