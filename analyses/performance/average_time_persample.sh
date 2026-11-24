#!/bin/bash
# Execution time per sample. PPR: the sum of the task realtimes in its Nextflow trace, without the structural-variant steps.
# GATK: the duration of the Nextflow run of the sample (from its report), typed in below. Output: total_time_results.csv.
# Usage: bash average_time_persample.sh [folder with PPR/<platform>/<sample>_trace.txt]

BASE_DIR=${1:-.}
cd $BASE_DIR

echo "=== CALCULATE TOTAL EXECUTION TIME PER SAMPLE ==="

# Function that calculates the total execution time (sum of realtime) - only for PPR
calculate_total_time_ppr() {
    local trace_file="$1"
    
    awk -F'\t' '
    BEGIN { total_seconds = 0 }
    
    NR > 1 {
        process_name = $4
        
        if (match(process_name, /^([A-Z_]+)/)) {
            process = substr(process_name, RSTART, RLENGTH)
        } else {
            process = process_name
        }
        
        # For PPR, exclude VG_CALL and VG_PACK (SV calling)
        if (process !~ /^(VG_CALL|VG_PACK|VG_DECONSTRUCT|VCFBUB_FILTER|VCFWAVE|PANGENIE_PREP|PANGENIE_INDEX|PANGENIE_GENOTYPE|PANGENIE_BIALLELIC|COMPRESS_INDEX_VCF_PANGENIE)$/) {
            time_str = $9  # realtime in column 9 for PPR
            
            if (time_str != "-" && time_str != "" && time_str != "realtime") {
                seconds = 0
                
                if (match(time_str, /([0-9]+)h/)) {
                    hours = substr(time_str, RSTART, RLENGTH-1)
                    seconds += hours * 3600
                }
                if (match(time_str, /([0-9]+)m/)) {
                    minutes = substr(time_str, RSTART, RLENGTH-1)
                    seconds += minutes * 60
                }
                if (match(time_str, /([0-9]+\.?[0-9]*)s/)) {
                    secs = substr(time_str, RSTART, RLENGTH-1)
                    seconds += int(secs + 0.5)
                }
                
                total_seconds += seconds
            }
        }
    }
    
    END { print total_seconds }
    ' "$trace_file"
}

# Function that converts seconds to a human-readable format
seconds_to_human() {
    local total_seconds=$1
    local hours=$((total_seconds / 3600))
    local minutes=$(((total_seconds % 3600) / 60))
    local seconds=$((total_seconds % 60))
    
    if [[ $hours -gt 0 ]]; then
        echo "${hours}h ${minutes}m ${seconds}s"
    elif [[ $minutes -gt 0 ]]; then
        echo "${minutes}m ${seconds}s"
    else
        echo "${seconds}s"
    fi
}

# The values for GATK Aviti from the Nextflow report of each run
declare -A gatk_aviti_hours=(
    ["1"]="15.17"   # 15h10m
    ["2"]="18.90"   # 18h54m
    ["3"]="14.85"   # 14h51m
    ["4"]="12.77"   # 12h46m
    ["5"]="10.78"   # 10h47m
    ["6"]="9.72"    # 9h43m
    ["7"]="8.73"    # 8h44m
    ["8"]="8.47"    # 8h28m
    ["9"]="13.27"   # 13h16m
    ["10"]="9.57"   # 9h34m
    ["11"]="10.65"  # 10h39m
    ["12"]="18.53"  # 18h32m
    ["13"]="10.40"  # 10h24m
    ["14"]="11.23"  # 11h14m
    ["15"]="11.23"  # 11h14m
    ["16"]="13.42"  # 13h25m
    ["17"]="12.20"  # 12h12m
    ["18"]="16.05"  # 16h03m
    ["19"]="9.52"   # 9h31m
    ["20"]="12.25"  # 12h15m
    ["21"]="19.73"  # 19h44m
    ["22"]="15.95"  # 15h57m
    ["24"]="7.23"   # 7h14m
)

# The values for GATK Illumina from the Nextflow report of each run (only the complete runs)
declare -A gatk_illumina_hours=(
    ["10"]="13.82"   # 13h49m
    ["12"]="12.65"  # 12h39m
    ["13"]="12.07"  # 12h04m
    ["14"]="14.97"  # 14h58m
    ["15"]="14.22"  # 14h13m
    ["16"]="11.07"  # 11h04m
    ["17"]="12.27"   # 12h16m
    ["18"]="14.42"   # 14h25m
    ["19"]="13.82"   # 13h49m
    ["20"]="12.45"  # 12h27m
    ["21"]="12.18"  # 12h11m
    ["22"]="8.22"   # 8h13m
    # Skip 23 with 1h19m (probably failed)
    ["24"]="12.32"  # 12h19m
    ["2"]="13.43"    # 13h26m
    ["3"]="11.77"    # 11h46m
    ["4"]="12.92"   # 12h55m
    ["5"]="13.75"   # 13h45m
    ["7"]="13.60"   # 13h36m
    ["8"]="38.25"   # 38h15m
    ["9"]="22.52"    # 22h31m
)

# Create the results file for the total time
results_file="total_time_results.csv"
echo "Pipeline,Sample,Total_Time_Hours" > $results_file

# Arrays that store the results
declare -a ppr_illumina_times
declare -a ppr_aviti_times
declare -a gatk_aviti_times
declare -a gatk_illumina_times

echo "Processing PPR Illumina samples..."
for trace_file in PPR/illumina/*_trace.txt; do
    if [[ -f "$trace_file" ]]; then
        sample_name=$(basename "$trace_file" _trace.txt)
        total_time=$(calculate_total_time_ppr "$trace_file")
        
        ppr_illumina_times+=($total_time)
        
        total_hours=$(echo "scale=2; $total_time / 3600" | bc)
        
        echo "  $sample_name: $(seconds_to_human $total_time)"
        echo "PPR_Illumina,$sample_name,$total_hours" >> $results_file
    fi
done

echo ""
echo "Processing PPR Aviti samples..."
for trace_file in PPR/aviti/*_trace.txt; do
    if [[ -f "$trace_file" ]]; then
        sample_name=$(basename "$trace_file" _trace.txt)
        total_time=$(calculate_total_time_ppr "$trace_file")
        
        ppr_aviti_times+=($total_time)
        
        total_hours=$(echo "scale=2; $total_time / 3600" | bc)
        
        echo "  $sample_name: $(seconds_to_human $total_time)"
        echo "PPR_Aviti,$sample_name,$total_hours" >> $results_file
    fi
done

echo ""
echo "Processing GATK Aviti samples (from the run reports)..."
for sample in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 24; do
    if [[ -n "${gatk_aviti_hours[$sample]}" ]]; then
        hours="${gatk_aviti_hours[$sample]}"
        seconds=$(echo "$hours * 3600" | bc | cut -d. -f1)
        gatk_aviti_times+=($seconds)
        
        echo "  $sample: ${hours}h"
        echo "GATK_Aviti,$sample,$hours" >> $results_file
    fi
done

echo ""
echo "Processing GATK Illumina samples (from the run reports)..."
for sample in 10 12 13 14 15 16 17 18 19 20 21 22 24 2 3 4 5 7 8 9; do
    if [[ -n "${gatk_illumina_hours[$sample]}" ]]; then
        hours="${gatk_illumina_hours[$sample]}"
        seconds=$(echo "$hours * 3600" | bc | cut -d. -f1)
        gatk_illumina_times+=($seconds)
        
        echo "  $sample: ${hours}h"
        echo "GATK_Illumina,$sample,$hours" >> $results_file
    fi
done

# Calculate the means
echo ""
echo "=== CALCULATE MEANS ==="

# PPR Illumina
if [[ ${#ppr_illumina_times[@]} -gt 0 ]]; then
    total_time=0
    for time in "${ppr_illumina_times[@]}"; do
        total_time=$((total_time + time))
    done
    avg_time_ppr_illumina=$((total_time / ${#ppr_illumina_times[@]}))
    avg_hours_ppr_illumina=$(echo "scale=2; $avg_time_ppr_illumina / 3600" | bc)
    
    echo "PPR Illumina - Mean time: $(seconds_to_human $avg_time_ppr_illumina)"
    echo "PPR_Illumina,AVERAGE,$avg_hours_ppr_illumina" >> $results_file
fi

# PPR Aviti
if [[ ${#ppr_aviti_times[@]} -gt 0 ]]; then
    total_time=0
    for time in "${ppr_aviti_times[@]}"; do
        total_time=$((total_time + time))
    done
    avg_time_ppr_aviti=$((total_time / ${#ppr_aviti_times[@]}))
    avg_hours_ppr_aviti=$(echo "scale=2; $avg_time_ppr_aviti / 3600" | bc)
    
    echo "PPR Aviti - Mean time: $(seconds_to_human $avg_time_ppr_aviti)"
    echo "PPR_Aviti,AVERAGE,$avg_hours_ppr_aviti" >> $results_file
fi

# GATK Aviti
if [[ ${#gatk_aviti_times[@]} -gt 0 ]]; then
    total_time=0
    for time in "${gatk_aviti_times[@]}"; do
        total_time=$((total_time + time))
    done
    avg_time_gatk_aviti=$((total_time / ${#gatk_aviti_times[@]}))
    avg_hours_gatk_aviti=$(echo "scale=2; $avg_time_gatk_aviti / 3600" | bc)
    
    echo "GATK Aviti - Mean time: $(seconds_to_human $avg_time_gatk_aviti)"
    echo "GATK_Aviti,AVERAGE,$avg_hours_gatk_aviti" >> $results_file
fi

# GATK Illumina
if [[ ${#gatk_illumina_times[@]} -gt 0 ]]; then
    total_time=0
    for time in "${gatk_illumina_times[@]}"; do
        total_time=$((total_time + time))
    done
    avg_time_gatk_illumina=$((total_time / ${#gatk_illumina_times[@]}))
    avg_hours_gatk_illumina=$(echo "scale=2; $avg_time_gatk_illumina / 3600" | bc)
    
    echo "GATK Illumina - Mean time: $(seconds_to_human $avg_time_gatk_illumina)"
    echo "GATK_Illumina,AVERAGE,$avg_hours_gatk_illumina" >> $results_file
fi

echo ""
echo "The results were saved to: $results_file"
