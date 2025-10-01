#!/bin/bash

# Configuration
BAM_DIR="/home/shsoudi/methylation/full_cohort_analysis/CAMDAC/bams_UPS/male/test/bams"
OUTPUT_DIR="/home/shsoudi/methylation/full_cohort_analysis/CAMDAC/bams_UPS/male/test/bams/CAMDAC_outs_UPS"
LOG_DIR="/home/shsoudi/methylation/full_cohort_analysis/CAMDAC/bams_UPS/male/test/bams/logs"  # New log directory
MERGED_DIR="${BAM_DIR}/merged_bams"
NORMAL_BAM="/home/shsoudi/methylation/full_cohort_analysis/CAMDAC/bams_UPS/proxy_normals/SRC345_N.sorted.bam"
PIPELINE_FILES="/home/shsoudi/softwares/CAMDAC_pipeline_files/pipeline_files_V1/"
LOG_FILE="${LOG_DIR}/merge_and_process.log"  # Logs go to log directory
COMPLETED_FILE="${LOG_DIR}/completed_samples.txt"  # Also in log directory

# Create directories
mkdir -p "$MERGED_DIR"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# R script template
R_SCRIPT_TEMPLATE=$(cat << 'EOF'
#!/usr/bin/env Rscript
library(CAMDAC)

args <- commandArgs(trailingOnly = TRUE)
patient_id <- args[1]
tumor_bam <- args[2]

tryCatch({
    CAMDAC::pipeline_tumor_normal(
        patient_id=patient_id,
        tumor_id="T",
        normal_id="N",
        tumor_bam=tumor_bam,
        normal_bam="NORMAL_BAM_PLACEHOLDER",
        sex="XY",
        path="OUTPUT_DIR_PLACEHOLDER",
        pipeline_files="PIPELINE_FILES_PLACEHOLDER",
        build="GRCH37",
        min_tumor = 1,
        min_normal = 1,
        mq = 0,
        n_cores = 30,
        paired_end = TRUE  
    )
    cat("SUCCESS:", patient_id, "\n")
}, error = function(e) {
    cat("ERROR:", patient_id, "-", conditionMessage(e), "\n")
    quit(save = "no", status = 1)
})
EOF
)

# Create log files
touch "$LOG_FILE"
touch "$COMPLETED_FILE"

# Function to log messages
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Function to check if sample is already completed
is_completed() {
    local patient_id="$1"
    if [[ -d "${OUTPUT_DIR}/${patient_id}" ]]; then
        return 0
    else
        return 1
    fi
}

# Function to detect all patient-timepoint groups and identify which need merging
detect_samples() {
    declare -A multi_region_groups
    declare -A single_region_samples
    
    # Find all BAM files
    for bam_file in "${BAM_DIR}"/*_pe_deduped.sorted.bam; do
        [ -e "$bam_file" ] || continue
        
        filename=$(basename "$bam_file")
        patient_id=$(basename "$bam_file" "_pe_deduped.sorted.bam")
        
        # Pattern: SRC127_T1_A_pe_deduped.sorted.bam → SRC127_T1
        if [[ "$filename" =~ ^(SRC[0-9]+_T[0-9]+)_[A-Z]_pe_deduped\.sorted\.bam$ ]]; then
            patient_timepoint="${BASH_REMATCH[1]}"
            
            # Count how many files for this patient_timepoint
            count=$(find "$BAM_DIR" -maxdepth 1 -name "${patient_timepoint}_*_pe_deduped.sorted.bam" -type f | wc -l)
            
            if [[ $count -gt 1 ]]; then
                multi_region_groups["$patient_timepoint"]=1
            else
                single_region_samples["$patient_id"]="$bam_file"
            fi
        fi
    done
    
    # Return results
    echo "MULTI_REGION:${!multi_region_groups[@]}"
    for sample_id in "${!single_region_samples[@]}"; do
        echo "SINGLE:$sample_id:${single_region_samples[$sample_id]}"
    done
}

# Function to merge BAM files for a patient-timepoint group
merge_patient_bams() {
    local patient_timepoint="$1"
    local merged_name="${patient_timepoint}_merged"
    
    # Skip if merged file already exists
    if [[ -f "${MERGED_DIR}/${merged_name}_pe_deduped.sorted.bam" ]]; then
        log_message "SKIPPING: Merged BAM already exists for $patient_timepoint"
        return 0
    fi
    
    # Find all BAM files for this patient_timepoint pattern
    local bam_files=($(find "$BAM_DIR" -maxdepth 1 -name "${patient_timepoint}_*_pe_deduped.sorted.bam" -type f | sort))
    
    if [[ ${#bam_files[@]} -gt 1 ]]; then
        log_message "MERGING: Found ${#bam_files[@]} BAM files for $patient_timepoint"
        log_message "FILES: ${bam_files[*]}"
        
        # Merge BAM files
        samtools merge -f "${MERGED_DIR}/${merged_name}_pe_deduped.sorted.bam" "${bam_files[@]}"
        
        if [[ $? -eq 0 ]]; then
            # Index the merged BAM
            samtools index "${MERGED_DIR}/${merged_name}_pe_deduped.sorted.bam"
            log_message "SUCCESS: Created ${MERGED_DIR}/${merged_name}_pe_deduped.sorted.bam"
            return 0
        else
            log_message "ERROR: BAM merge failed for $patient_timepoint"
            return 1
        fi
    else
        log_message "SKIPPING: Only 1 BAM file found for $patient_timepoint"
        return 0
    fi
}

# Function to process a single sample with CAMDAC
process_with_camdac() {
    local bam_file="$1"
    local patient_id="$2"
    
    # Skip if output already exists
    if is_completed "$patient_id"; then
        log_message "SKIPPING CAMDAC: Output already exists for $patient_id"
        return 0
    fi
    
    log_message "PROCESSING CAMDAC: $patient_id"
    
    # Create temporary R script
    local temp_r_script=$(mktemp)
    echo "$R_SCRIPT_TEMPLATE" | \
        sed "s|NORMAL_BAM_PLACEHOLDER|${NORMAL_BAM}|g" | \
        sed "s|OUTPUT_DIR_PLACEHOLDER|${OUTPUT_DIR}|g" | \
        sed "s|PIPELINE_FILES_PLACEHOLDER|${PIPELINE_FILES}|g" > "$temp_r_script"
    
    # Execute R script - save individual logs to log directory
    local output_file="${LOG_DIR}/${patient_id}_camdac.log"
    if Rscript "$temp_r_script" "$patient_id" "$bam_file" 2>&1 | tee -a "$output_file"; then
        log_message "SUCCESS CAMDAC: $patient_id"
    else
        log_message "FAILED CAMDAC: $patient_id"
    fi
    
    rm -f "$temp_r_script"
}

# Main execution
log_message "=== Starting automated merge and CAMDAC processing ==="
log_message "BAM Directory: $BAM_DIR"
log_message "Output Directory: $OUTPUT_DIR"
log_message "Log Directory: $LOG_DIR"
log_message "Merged BAMs Directory: $MERGED_DIR"

# Step 1: Detect which samples need merging and which don't
log_message "Step 1: Detecting sample types..."
detection_results=$(detect_samples)

multi_region_groups=()
declare -A single_region_samples

while IFS= read -r line; do
    if [[ "$line" == MULTI_REGION:* ]]; then
        multi_region_groups=(${line#MULTI_REGION:})
    elif [[ "$line" == SINGLE:* ]]; then
        IFS=':' read -r prefix sample_id bam_file <<< "$line"
        single_region_samples["$sample_id"]="$bam_file"
    fi
done <<< "$detection_results"

log_message "Multi-region groups (will be merged): ${#multi_region_groups[@]} - ${multi_region_groups[*]}"
log_message "Single-region samples (will be processed directly): ${#single_region_samples[@]} - ${!single_region_samples[*]}"

# Step 2: Merge multi-region samples
log_message "Step 2: Merging multi-region samples..."
for group in "${multi_region_groups[@]}"; do
    merge_patient_bams "$group"
done

# Step 3: Run CAMDAC on merged samples
log_message "Step 3: Running CAMDAC on merged samples..."
for merged_bam in "${MERGED_DIR}"/*_merged_pe_deduped.sorted.bam; do
    [ -e "$merged_bam" ] || continue
    
    patient_id=$(basename "$merged_bam" "_merged_pe_deduped.sorted.bam")
    log_message "Processing merged sample: $patient_id"
    process_with_camdac "$merged_bam" "$patient_id"
done

# Step 4: Run CAMDAC on single-region samples
log_message "Step 4: Running CAMDAC on single-region samples..."
for sample_id in "${!single_region_samples[@]}"; do
    bam_file="${single_region_samples[$sample_id]}"
    if [[ -f "$bam_file" ]]; then
        log_message "Processing single-region sample: $sample_id"
        process_with_camdac "$bam_file" "$sample_id"
    else
        log_message "ERROR: BAM file not found for $sample_id: $bam_file"
    fi
done

log_message "=== Merge and CAMDAC processing completed ==="
log_message "Multi-region samples processed: ${#multi_region_groups[@]}"
log_message "Single-region samples processed: ${#single_region_samples[@]}"
log_message "Main log file: $LOG_FILE"
log_message "Individual sample logs: ${LOG_DIR}/*_camdac.log"
log_message "Completed samples list: $COMPLETED_FILE"
log_message "Merged BAMs are in: $MERGED_DIR"
log_message "CAMDAC output folders are in: $OUTPUT_DIR"
