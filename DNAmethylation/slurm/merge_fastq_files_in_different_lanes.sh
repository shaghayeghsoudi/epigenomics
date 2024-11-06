#!/bin/bash -l
# ---------------------------------------------------------------------
# SLURM script for demultiplexed fumi_tools
# ---------------------------------------------------------------------
#SBATCH --job-name==fumi_tools
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=14G
#SBATCH --tasks=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shsoudi
#SBATCH --partition=emoding
#SBATCH --export=ALL




# Define the folder containing your fastq.gz files
FOLDER="/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/NovaSeqX15-22_processing_upstream"

# Move into the folder
cd "$FOLDER" || exit 1  # Exit if the folder doesn't exist

# Identify prefixes with repeated files and concatenate them
for prefix in $(ls *_R1.fastq.gz | sed -E 's/_L00[0-9]+_R1.fastq.gz//' | sort | uniq -d); do
    # Collect all R1 files with the repeated prefix
    R1_files=(${prefix}_L00*_R1.fastq.gz)
    
    # Define the output filename for the merged file
    R1_merged="${prefix}_R1_merged.fastq.gz"

    # Concatenate R1 files
    cat "${R1_files[@]}" > "$R1_merged"
    echo "Merged R1 files for $prefix into $R1_merged"

    rm "${R1_files[@]}"
    echo "Removed original R1 files for $prefix"
done


for prefix in $(ls *_R2.fastq.gz | sed -E 's/_L00[0-9]+_R2.fastq.gz//' | sort | uniq -d); do
    # Collect all R1 files with the repeated prefix
    R2_files=(${prefix}_L00*_R2.fastq.gz)
    
    # Define the output filename for the merged file
    R2_merged="${prefix}_R2_merged.fastq.gz"

    # Concatenate R1 files
    cat "${R2_files[@]}" > "$R2_merged"
    echo "Merged R2 files for $prefix into $R2_merged"


    rm "${R2_files[@]}"
    echo "Removed original R2 files for $prefix"


done

