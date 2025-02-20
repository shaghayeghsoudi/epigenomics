#!/bin/bash -l
# ---------------------------------------------------------------------
# SLURM script for spike_processing
# ---------------------------------------------------------------------
#SBATCH --job-name==spike_processing
#SBATCH --cpus-per-task=8
#SBATCH --array=1-92
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shsoudi
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


## Directories:
# Locate the input data
ROOT_DIR=/oak/stanford/groups/emoding/analysis/shaghayegh/methylation/NovaSeqX15-22_processing_upstream

### REF directories
REF_DIR_Meth=/oak/stanford/groups/emoding/analysis/shaghayegh/resources/GRCh37/Bisulfite_Genome_spike_control/RRBS_methylated_control
REF_DIR_Unmeth=/oak/stanford/groups/emoding/analysis/shaghayegh/resources/GRCh37/Bisulfite_Genome_spike_control/RRBS_unmethylated_control

CONFIG=${ROOT_DIR}/sample_ID_config.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
BASE_NAME=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $CONFIG)
SAMPLE_R1=${ROOT_DIR}/trimmed_fastq/${BASE_NAME}_R1_val_1.fq.gz    
SAMPLE_R2=${ROOT_DIR}/trimmed_fastq/${BASE_NAME}_R2_val_2.fq.gz    


##########################################################
####### Step 1: Activate Bismark and do Alignment ########
##########################################################

echo "Activating bismark environment..."
source ~/anaconda3/etc/profile.d/conda.sh
source activate /home/users/shsoudi/emoding/envs/bismark
echo "Starting bismark alignment..."


#OUTPUT_DIR=${ROOT_DIR}/NovaSeqX10_P1/NovaSeqX10_P1_alignmnet_bismark_trimmed
OUTPUT_DIR_BIS=${ROOT_DIR}/alignmnet_bismark_scored_trimmed_spiked
mkdir -p ${OUTPUT_DIR_BIS}

#TEMP_DIR=${ROOT_DIR}/NovaSeqX10_P1/NovaSeqX10_P1_alignmnet_bismark_temp
TEMP_DIR=${ROOT_DIR}/alignmnet_bismark_temp_scored_trimmed_spiked   
mkdir -p ${TEMP_DIR}

#bismark --pbat -1 ${SAMPLE_R1} -2 ${SAMPLE_R2} --bowtie2 --bam --temp_dir ${TEMP_DIR} --genome ${REF_DIR} -o ${OUTPUT_DIR} --un

### without spikein 
#bismark --pbat --score_min L,0,-0.6 -1 ${SAMPLE_R1} -2 ${SAMPLE_R2} --bowtie2 --bam --temp_dir ${TEMP_DIR} --genome ${REF_DIR} -o ${OUTPUT_DIR} --un

### with spike control
bismark -q --pbat --score_min L,0,-0.6 --prefix Meth_ctrl ${REF_DIR_Meth} -1 ${SAMPLE_R1} -2 ${SAMPLE_R2} --bowtie2 --bam --temp_dir ${TEMP_DIR} -o ${OUTPUT_DIR_BIS} --un
bismark -q --pbat --score_min L,0,-0.6 --prefix unmeth_ctrl ${REF_DIR_Unmeth} -1 ${SAMPLE_R1} -2 ${SAMPLE_R2} --bowtie2 --bam --temp_dir ${TEMP_DIR} -o ${OUTPUT_DIR_BIS} --un

echo "bismark alignment completed."
conda deactivate  # Deactivate BWA environment


##########################################
### Step 2: indexing and deduplicating ###
##########################################

echo "Activating Samtools environment..."
source ~/anaconda3/etc/profile.d/conda.sh
source activate /home/users/shsoudi/emoding/envs/samtools112
echo "Running samtools indexing and sorting..."



# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
SAMPLE_BAM_Meth=${OUTPUT_DIR_BIS}/Meth_ctrl.${BASE_NAME}_R1_val_1_bismark_bt2_pe.bam 
SAMPLE_BAM_UnMeth=${OUTPUT_DIR_BIS}/unmeth_ctrl.${BASE_NAME}_R1_val_1_bismark_bt2_pe.bam 


OUTPUT_DIR_DEDUP=${OUTPUT_DIR_BIS}/dedupped_bams_fumitools
mkdir -p ${OUTPUT_DIR_DEDUP}

### methylated bam
samtools sort ${SAMPLE_BAM_Meth} -o ${OUTPUT_DIR_DEDUP}/Meth_ctrl.${BASE_NAME}_initial_sorted_undedupped.bam 
### unmethylated bam
samtools sort ${SAMPLE_BAM_UnMeth} -o ${OUTPUT_DIR_DEDUP}/UnMeth_ctrl.${BASE_NAME}_initial_sorted_undedupped.bam 

echo "samtools index completed."
echo "Running fumi_tools..."

fumi_tools dedup -i ${OUTPUT_DIR_DEDUP}/Meth_ctrl.${BASE_NAME}_initial_sorted_undedupped.bam  -o ${OUTPUT_DIR_DEDUP}/Meth_ctrl.${BASE_NAME}_deduped.bam --paired
fumi_tools dedup -i ${OUTPUT_DIR_DEDUP}/UnMeth_ctrl.${BASE_NAME}_initial_sorted_undedupped.bam  -o ${OUTPUT_DIR_DEDUP}/UnMeth_ctrl.${BASE_NAME}_deduped.bam --paired


echo "fumi_tools dedup completed."
conda deactivate  # Deactivate BWA environment



####################################################################
#### Step 3: Activate Bismark conda and do methylation calling #####
####################################################################

echo "Activating bismark environment..."
source ~/anaconda3/etc/profile.d/conda.sh
source activate /home/users/shsoudi/emoding/envs/bismark
echo "Running bismark methylation calling..."


BAM_Meth_DEDUP=${OUTPUT_DIR_DEDUP}/Meth_ctrl.${BASE_NAME}_deduped.bam
BAM_UnMeth_DEDUP=${OUTPUT_DIR_DEDUP}/UnMeth_ctrl.${BASE_NAME}_deduped.bam


OUTPUT_DIR_CALL=${ROOT_DIR}/methylation_calls_cli_spiked
mkdir -p ${OUTPUT_DIR}


bismark_methylation_extractor -p --gzip --bedGraph --CX ${BAM_Meth_DEDUP} -o ${OUTPUT_DIR_CALL}
bismark_methylation_extractor -p --gzip --bedGraph --CX ${BAM_UnMeth_DEDUP} -o ${OUTPUT_DIR_CALL}

#### Bismark2report
#BAM_REPORT=${ROOT_DIR}/NovaSeqX10_P1/NovaSeqX10_P1_alignmnet_bismark_scored0.6_trimmed/${BASE_NAME}_L006_R1_val_1_bismark_bt2_PE_report.txt
#mbias_report=${OUTPUT_DIR}/${BASE_NAME}_deduped.M-bias.txt
#splitting_report=${OUTPUT_DIR}/${BASE_NAME}_deduped_splitting_report.txt

#bismark2report --dir ${OUTPUT_DIR} --alignment_report ${BAM_REPORT} --splitting_report ${splitting_report} --mbias_report ${mbias_report} 

echo "bismark methylation calling completed."
conda deactivate  # Deactivate BWA environment
