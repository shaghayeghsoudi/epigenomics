################################
################################

#!/bin/bash -l
# ---------------------------------------------------------------------
# SLURM script for demultiplexed
# ---------------------------------------------------------------------
#SBATCH --job-name==bbmap
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=50G
#SBATCH --tasks=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shsoudi
#SBATCH --partition=emoding
#SBATCH --export=ALL
#module purge
# Activate Anaconda work environment for OpenDrift

module load java/17.0.4


## Directories:
# Locate the input data
SCRIPT_DIR=/oak/stanford/groups/emoding/analysis/shaghayegh/softwares/bbmap
ROOT_DIR=/oak/stanford/groups/emoding/sequencing/methylation_RRBS/p1/usftp21.novogene.com/01.RawData
FASQ_DIR_RAW=${ROOT_DIR}/P1


DEMULTIPLEX_DIR=${ROOT_DIR}/P1/demultiplexed



bash ${SCRIPT_DIR}/demuxbyname.sh -Xmx10g 
in1=${FASQ_DIR_RAW}/P1_CKDL230029542-1A_22CGCNLT3_L3_1.fq.gz 
in2=${FASQ_DIR_RAW}/P1_CKDL230029542-1A_22CGCNLT3_L3_2.fq.gz 
out1=${DEMULTIPLEX_DIR}/out_%_R1.fq.gz 
out2=${DEMULTIPLEX_DIR}/out_%_R2.fq.gz names=${FASQ_DIR_RAW}/names.txt 
threads=$SLURM_CPUS_PER_TASK


