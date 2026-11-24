#!/bin/bash
# CPU minutes per sample (realtime x %cpu summed over the tasks of the Nextflow trace), without the QC steps of GATK and the
# structural-variant steps of PPR. Output: cpu_minutes_results.csv.
# Usage: bash cpumin.sh [folder with PPR/<platform>/ and GATK/<platform>/individual/ traces]

BASE_DIR=${1:-.}
cd $BASE_DIR

echo "=== CALCULATE CPU MINUTES PER SAMPLE ==="

# Function that calculates the CPU minutes from a trace file
calculate_cpu_minutes() {
    local trace_file="$1"
    local pipeline_type="$2"
    
    # Find the columns for the pipeline type
    if [[ "$pipeline_type" == "gatk_aviti" ]]; then
        process_col=4
        cpu_col=22
        realtime_col=20
    else
        process_col=4
        cpu_col=10
        realtime_col=9
    fi
    
    awk -F'\t' -v proc_col="$process_col" -v cpu_col="$cpu_col" -v rt_col="$realtime_col" '
    BEGIN { total_cpu_minutes = 0 }
    
    # Skip the header and the QC processes
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
            cpu_percent = $cpu_col
            realtime = $rt_col
            
            if (cpu_percent != "-" && cpu_percent != "" && cpu_percent != "%cpu" && 
                realtime != "-" && realtime != "" && realtime != "realtime") {
                
                # Convert realtime to minutes
                realtime_minutes = 0
                if (match(realtime, /([0-9]+)h/)) {
                    hours = substr(realtime, RSTART, RLENGTH-1)
                    realtime_minutes += hours * 60
                }
                if (match(realtime, /([0-9]+)m/)) {
                    minutes = substr(realtime, RSTART, RLENGTH-1)
                    realtime_minutes += minutes
                }
                if (match(realtime, /([0-9]+\.?[0-9]*)s/)) {
                    seconds = substr(realtime, RSTART, RLENGTH-1)
                    realtime_minutes += seconds / 60
                }
                
                # Convert the CPU percent (for example "5933.1%") to a number
                gsub(/%/, "", cpu_percent)
                cpu_decimal = cpu_percent / 100
                
                # Calculate CPU minutes = realtime_minutes * cpu_utilization
                cpu_minutes = realtime_minutes * cpu_decimal
                total_cpu_minutes += cpu_minutes
            }
        }
    }
    
    END { print total_cpu_minutes }
    ' "$trace_file"
}

# Create the results file for the CPU minutes
results_file="cpu_minutes_results.csv"
echo "Pipeline,Sample,CPU_Minutes" > $results_file

# Arrays that store the results
declare -a ppr_illumina_cpu
declare -a ppr_aviti_cpu
declare -a gatk_aviti_cpu

echo "Processing PPR Illumina samples..."
for trace_file in PPR/illumina/*_trace.txt; do
    if [[ -f "$trace_file" ]]; then
        sample_name=$(basename "$trace_file" _trace.txt)
        cpu_minutes=$(calculate_cpu_minutes "$trace_file" "ppr")
        
        ppr_illumina_cpu+=($cpu_minutes)
        
        cpu_min_rounded=$(echo "scale=1; $cpu_minutes" | bc)
        
        echo "  $sample_name: ${cpu_min_rounded} CPU minutes"
        echo "PPR_Illumina,$sample_name,$cpu_min_rounded" >> $results_file
    fi
done

echo ""
echo "Processing PPR Aviti samples..."
for trace_file in PPR/aviti/*_trace.txt; do
    if [[ -f "$trace_file" ]]; then
        sample_name=$(basename "$trace_file" _trace.txt)
        cpu_minutes=$(calculate_cpu_minutes "$trace_file" "ppr")
        
        ppr_aviti_cpu+=($cpu_minutes)
        
        cpu_min_rounded=$(echo "scale=1; $cpu_minutes" | bc)
        
        echo "  $sample_name: ${cpu_min_rounded} CPU minutes"
        echo "PPR_Aviti,$sample_name,$cpu_min_rounded" >> $results_file
    fi
done

echo ""
echo "Processing GATK Aviti samples..."
for trace_file in GATK/aviti/individual/*_trace.txt; do
    if [[ -f "$trace_file" ]]; then
        sample_name=$(basename "$trace_file" _trace.txt)
        cpu_minutes=$(calculate_cpu_minutes "$trace_file" "gatk_aviti")
        
        gatk_aviti_cpu+=($cpu_minutes)
        
        cpu_min_rounded=$(echo "scale=1; $cpu_minutes" | bc)
        
        echo "  $sample_name: ${cpu_min_rounded} CPU minutes"
        echo "GATK_Aviti,$sample_name,$cpu_min_rounded" >> $results_file
    fi
done

# Calculate the means
echo ""
echo "=== CALCULATE MEAN CPU MINUTES ==="

# PPR Illumina
if [[ ${#ppr_illumina_cpu[@]} -gt 0 ]]; then
    total_cpu=0
    for cpu in "${ppr_illumina_cpu[@]}"; do
        total_cpu=$(echo "$total_cpu + $cpu" | bc)
    done
    avg_cpu_ppr_illumina=$(echo "scale=1; $total_cpu / ${#ppr_illumina_cpu[@]}" | bc)
    
    echo "PPR Illumina - Mean CPU: ${avg_cpu_ppr_illumina} minutes"
    echo "PPR_Illumina,AVERAGE,$avg_cpu_ppr_illumina" >> $results_file
fi

# PPR Aviti
if [[ ${#ppr_aviti_cpu[@]} -gt 0 ]]; then
    total_cpu=0
    for cpu in "${ppr_aviti_cpu[@]}"; do
        total_cpu=$(echo "$total_cpu + $cpu" | bc)
    done
    avg_cpu_ppr_aviti=$(echo "scale=1; $total_cpu / ${#ppr_aviti_cpu[@]}" | bc)
    
    echo "PPR Aviti - Mean CPU: ${avg_cpu_ppr_aviti} minutes"
    echo "PPR_Aviti,AVERAGE,$avg_cpu_ppr_aviti" >> $results_file
fi

# GATK Aviti
if [[ ${#gatk_aviti_cpu[@]} -gt 0 ]]; then
    total_cpu=0
    for cpu in "${gatk_aviti_cpu[@]}"; do
        total_cpu=$(echo "$total_cpu + $cpu" | bc)
    done
    avg_cpu_gatk_aviti=$(echo "scale=1; $total_cpu / ${#gatk_aviti_cpu[@]}" | bc)
    
    echo "GATK Aviti - Mean CPU: ${avg_cpu_gatk_aviti} minutes"
    echo "GATK_Aviti,AVERAGE,$avg_cpu_gatk_aviti" >> $results_file
fi

echo ""
echo "The results were saved to: $results_file"

echo ""
echo "Processing GATK Illumina samples..."
for trace_file in GATK/illumina/individual/*_trace.txt; do
    if [[ -f "$trace_file" ]]; then
        sample_name=$(basename "$trace_file" _trace.txt)
        # Skip the problem sample
        if [[ "$sample_name" == "23" ]]; then
            continue
        fi

        cpu_minutes=$(calculate_cpu_minutes "$trace_file" "gatk_illumina")

        gatk_illumina_cpu+=($cpu_minutes)

        cpu_min_rounded=$(echo "scale=1; $cpu_minutes" | bc)

        echo "  $sample_name: ${cpu_min_rounded} CPU minutes"
        echo "GATK_Illumina,$sample_name,$cpu_min_rounded" >> $results_file
    fi
done

# GATK Illumina
if [[ ${#gatk_illumina_cpu[@]} -gt 0 ]]; then
    total_cpu=0
    for cpu in "${gatk_illumina_cpu[@]}"; do
        total_cpu=$(echo "$total_cpu + $cpu" | bc)
    done
    avg_cpu_gatk_illumina=$(echo "scale=1; $total_cpu / ${#gatk_illumina_cpu[@]}" | bc)

    echo "GATK Illumina - Mean CPU: ${avg_cpu_gatk_illumina} minutes"
    echo "GATK_Illumina,AVERAGE,$avg_cpu_gatk_illumina" >> $results_file
fi


