#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=3G
#SBATCH --job-name=01-visium_image
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -o ../../processed-data/nda-submission/logs/01-visium_image.txt
#SBATCH -e ../../processed-data/nda-submission/logs/01-visium_image.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.3.x

## List current modules for reproducibility
module list

Rscript 01-visium_image.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/
