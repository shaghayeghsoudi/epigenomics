#!/bin/bash

### convert DSS to coverage
# Set input and output directories
INPUT_DIR="/home/shsoudi/methylation/full_cohort_analysis/DSS_input_filt_5reads_merged_per_subject/processed_file_merged_by_addingupp"  # Change this to your DSS folder path
OUTPUT_DIR="/home/shsoudi/methylation/full_cohort_analysis/coverage_files_merged_regions"    # Change this to your output folder path

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all DSS files in the input directory
for file in "$INPUT_DIR"/*.txt; do
    # Define output file name
    filename=$(basename "$file")
    output_file="$OUTPUT_DIR/${filename}_bismark.cov"

    # Process file and convert format
    awk '{
        end = $2;  # End position is same as Start
        if (($3 + $4) > 0) {
            meth_percentage = ($3 / ($3 + $4)) * 100;
        } else {
            meth_percentage = 0;
        }
        printf "%s\t%s\t%s\t%.2f\t%s\t%s\n", $1, $2, end, meth_percentage, $3, $4;
    }' "$file" > "$output_file"

    echo "Converted: $file -> $output_file"
done

echo "Conversion completed for all files!"