#!/bin/bash
#
#SBATCH --job-name=wget_fastq
#
#SBATCH --time=02-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shsoudi
#SBATCH --partition=emoding

mkdir -p /oak/stanford/groups/emoding/sequencing/methylation_RRBS/p1
cd /oak/stanford/groups/emoding/sequencing/methylation_RRBS/p1
wget -r -c ftp://X202SC23090732-Z01-F001:yrcahf2y@usftp21.novogene.com:21/

