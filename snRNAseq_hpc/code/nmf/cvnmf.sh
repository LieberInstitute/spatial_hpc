#!/bin/bash
#SBATCH --job-name=cvnmf
#SBATCH --output=logs/cvnmf.txt
#SBATCH --error=logs/cvnmf.txt
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=erik.nelson116@gmail.com # Please replace with the appropriate email address
echo "**** Job starts ****"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.3

## List current modules for reproducibility
module list

## Edit with your job command
Rscript cvnmf.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
