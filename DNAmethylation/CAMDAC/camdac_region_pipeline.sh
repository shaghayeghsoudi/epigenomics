#!/bin/bash

# Configuration
BAM_DIR="/home/shsoudi/methylation/full_cohort_analysis/CAMDAC/bams_UPS/male"
OUTPUT_DIR="/home/shsoudi/methylation/full_cohort_analysis/CAMDAC/CAMDAC_outs_UPS"
NORMAL_BAM="/home/shsoudi/methylation/full_cohort_analysis/CAMDAC/bams_UPS/proxy_normals/SRC345_N.sorted.bam"
PIPELINE_FILES="/home/shsoudi/softwares/CAMDAC_pipeline_files/pipeline_files_V1/"
LOG_FILE="${BAM_DIR}/camdac_batch_processing.log"
COMPLETED_FILE="${BAM_DIR}/completed_samples.txt"

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
        n_cores = 20,
        paired_end = TRUE  
    )
    cat("SUCCESS:", patient_id, "\n")
}, error = function(e) {
    cat("ERROR:", patient_id, "-", conditionMessage(e), "\n")
    quit(save = "no", status = 1)
})
EOF
)

# Create output directory and log files
mkdir -p "$OUTPUT_DIR"
touch "$LOG_FILE"
touch "$COMPLETED_FILE"

# Function to log messages
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Function to check if sample is already completed
is_completed() {
    local patient_id="$1"
    # Check if the output folder exists in the NEW location
    if [[ -d "${OUTPUT_DIR}/${patient_id}" ]]; then
        return 0
    else
        return 1
    fi
}

# Function to mark sample as completed
mark_completed() {
    echo "$1" >> "$COMPLETED_FILE"
}

# Function to process a single sample
process_sample() {
    local bam_file="$1"
    local patient_id=$(basename "$bam_file" "_pe_deduped.sorted.bam")
    
    # Skip if already completed (CAMDAC output folder exists in NEW location)
    if is_completed "$patient_id"; then
        log_message "SKIPPING: $patient_id - output folder already exists"
        return 0
    fi
    
    log_message "STARTING: Processing $patient_id"
    log_message "CAMDAC will create: ${OUTPUT_DIR}/${patient_id}/"
    
    # Create temporary R script
    local temp_r_script=$(mktemp)
    echo "$R_SCRIPT_TEMPLATE" | \
        sed "s|NORMAL_BAM_PLACEHOLDER|${NORMAL_BAM}|g" | \
        sed "s|OUTPUT_DIR_PLACEHOLDER|${OUTPUT_DIR}|g" | \
        sed "s|PIPELINE_FILES_PLACEHOLDER|${PIPELINE_FILES}|g" > "$temp_r_script"  # ← REMOVED THE STRAY BACKSLASH
    
    # Execute R script and capture output
    local output_file="${BAM_DIR}/${patient_id}_camdac_output.log"
    
    if Rscript "$temp_r_script" "$patient_id" "$bam_file" 2>&1 | tee -a "$output_file" | tee -a "$LOG_FILE"; then
        log_message "COMPLETED: Successfully processed $patient_id"
        mark_completed "$patient_id"
    else
        log_message "FAILED: Error processing $patient_id - check ${output_file} for details"
    fi
    
    # Clean up temporary script
    rm -f "$temp_r_script"
}

# Main processing loop
log_message "=== Starting CAMDAC batch processing ==="
log_message "BAM Directory: $BAM_DIR"
log_message "Output Directory: $OUTPUT_DIR"

# Process all BAM files
for bam_file in "${BAM_DIR}"/*_pe_deduped.sorted.bam; do
    [ -e "$bam_file" ] || continue
    process_sample "$bam_file"
done

log_message "=== Batch processing completed ==="

# Summary
log_message "=== Processing Summary ==="
log_message "Total samples completed: $(wc -l < "$COMPLETED_FILE" 2>/dev/null || echo 0)"
