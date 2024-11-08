python
Copy code
# Define global variables
FOLDER = "path/to/your/folder"  # Folder containing raw fastq.gz files

# List unique prefixes based on repeated files before "L00" for both R1 and R2
PREFIXES_R1 = [prefix for prefix in
               shell(f"cd {FOLDER} && ls *_R1.fastq.gz | sed -E 's/_L00[0-9]+_R1.fastq.gz//' | sort | uniq -d",
                     read=True).strip().split('\n')]

PREFIXES_R2 = [prefix for prefix in
               shell(f"cd {FOLDER} && ls *_R2.fastq.gz | sed -E 's/_L00[0-9]+_R2.fastq.gz//' | sort | uniq -d",
                     read=True).strip().split('\n')]

# Define default target rule
rule all:
    input:
        expand("{prefix}_R1_merged.fastq.gz", prefix=PREFIXES_R1),
        expand("{prefix}_R2_merged.fastq.gz", prefix=PREFIXES_R2)

# Rule to merge R1 files with repeated prefixes
rule merge_repeated_R1:
    input:
        lambda wildcards: expand(FOLDER + "/{prefix}_L00*_R1.fastq.gz", prefix=wildcards.prefix)
    output:
        "{prefix}_R1_merged.fastq.gz"
    shell:
        """
        cat {input} > {output}
        echo "Merged R1 files for {wildcards.prefix} into {output}"
        """

# Rule to merge R2 files with repeated prefixes
rule merge_repeated_R2:
    input:
        lambda wildcards: expand(FOLDER + "/{prefix}_L00*_R2.fastq.gz", prefix=wildcards.prefix)
    output:
        "{prefix}_R2_merged.fastq.gz"
    shell:
        """
        cat {input} > {output}
        echo "Merged R2 files for {wildcards.prefix} into {output}"
        """





##########
### remove repeated files 


import glob
import os

# Define input folder and files
FOLDER = "path/to/your/folder"
R1_files = glob.glob(os.path.join(FOLDER, "*_R1.fastq.gz"))
R2_files = glob.glob(os.path.join(FOLDER, "*_R2.fastq.gz"))

# Find repeated prefixes for R1
R1_prefixes = [os.path.basename(f).split("_L00")[0] for f in R1_files]
repeated_R1_prefixes = [prefix for prefix in set(R1_prefixes) if R1_prefixes.count(prefix) > 1]

# Find repeated prefixes for R2
R2_prefixes = [os.path.basename(f).split("_L00")[0] for f in R2_files]
repeated_R2_prefixes = [prefix for prefix in set(R2_prefixes) if R2_prefixes.count(prefix) > 1]

# Define targets for merged files
R1_targets = [os.path.join(FOLDER, f"{prefix}_R1_merged.fastq.gz") for prefix in repeated_R1_prefixes]
R2_targets = [os.path.join(FOLDER, f"{prefix}_R2_merged.fastq.gz") for prefix in repeated_R2_prefixes]

# Define main rule to merge and clean up
rule all:
    input:
        R1_targets + R2_targets

# Rule to concatenate and clean up R1 files for each repeated prefix
rule merge_and_cleanup_R1_files:
    input:
        lambda wildcards: sorted(glob.glob(os.path.join(FOLDER, f"{wildcards.prefix}_L00*_R1.fastq.gz")))
    output:
        merged=os.path.join(FOLDER, "{prefix}_R1_merged.fastq.gz")
    shell:
        """
        cat {input} > {output.merged}
        echo "Merged R1 files for {wildcards.prefix} into {output.merged}"
        rm {input}
        echo "Removed original R1 files for {wildcards.prefix}"
        """

# Rule to concatenate and clean up R2 files for each repeated prefix
rule merge_and_cleanup_R2_files:
    input:
        lambda wildcards: sorted(glob.glob(os.path.join(FOLDER, f"{wildcards.prefix}_L00*_R2.fastq.gz")))
    output:
        merged=os.path.join(FOLDER, "{prefix}_R2_merged.fastq.gz")
    shell:
        """
        cat {input} > {output.merged}
        echo "Merged R2 files for {wildcards.prefix} into {output.merged}"
        rm {input}
        echo "Removed original R2 files for {wildcards.prefix}"
        """


