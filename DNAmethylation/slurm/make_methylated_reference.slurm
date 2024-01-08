#!/bin/bash -l
# ---------------------------------------------------------------------
# SLURM script for trim_glore
# ---------------------------------------------------------------------
#SBATCH --job-name==trim_glore
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=20G
#SBATCH --tasks=1
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shsoudi
#SBATCH --partition=emoding
#SBATCH --export=ALL
#module purge
# Activate Anaconda work environment for OpenDrift

source ~/anaconda3/etc/profile.d/conda.sh
source activate /home/users/shsoudi/emoding/envs/bismark 



/home/users/shsoudi/emoding/envs/bismark/bin/bismark_genome_preparation --path_to_aligner 
/home/users/shsoudi/emoding/envs/bismark/bin/ --verbose 
/oak/stanford/groups/emoding/analysis/shaghayegh/resources/reference_GRCh37/

