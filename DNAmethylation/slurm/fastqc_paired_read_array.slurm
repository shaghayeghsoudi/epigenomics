#!/bin/bash -l
# ---------------------------------------------------------------------
# SLURM script for fastqc
# ---------------------------------------------------------------------
#SBATCH --job-name==fastqc
#SBATCH --cpus-per-task=2
#SBATCH --array=1-3
#SBATCH --mem-per-cpu=40G
#SBATCH --tasks=1
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shsoudi
#SBATCH --partition=emoding
#SBATCH --export=ALL
#module purge
# Activate Anaconda work environment for OpenDrift


# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------
echo ""
echo "Job Array ID / Job ID: $SLURM_ARRAY_JOB_ID / $SLURM_JOB_ID"
echo "This is job $SLURM_ARRAY_TASK_ID out of $SLURM_ARRAY_TASK_COUNT jobs."
echo ""

## Locate the fastqc
source ~/anaconda3/etc/profile.d/conda.sh
source activate /home/users/shsoudi/emoding/envs/fastqc


## Directories:
# Locate the input data
ROOT_DIR=/oak/stanford/groups/emoding/sequencing/methylation_RRBS/p1/usftp21.novogene.com/01.RawData/P1

# Specify the path to the config file
CONFIG=${ROOT_DIR}/demultiplexed/sample_ID_config.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
BASE_NAME=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIG)
SAMPLE=${ROOT_DIR}/demultiplexed/${BASE_NAME}


OUTPUT_DIR=${ROOT_DIR}/demultiplexed/fastqc_reports



fastqc ${SAMPLE}_R1.fq.gz ${SAMPLE}_R2.fq.gz -o ${OUTPUT_DIR}


# ---------------------------------------------------------------------
echo "Job finished with exit code $? at: `date`"
# ---------------------------------------------------------------------

