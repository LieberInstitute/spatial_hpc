#!/bin/bash
#SBATCH --job-name=nnSVG_PRECAST_batch_k5-20
#SBATCH --output=logs/nnsvg_PRECAST_batch_k.%a.txt
#SBATCH --error=logs/nnsvg_PRECAST_batch_k.%a.txt
#SBATCH --mail-type=END
#SBATCH --array=5-20%5
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=enelso40@jhmi.edu # Please replace with the appropriate email address
echo "**** Job starts ****"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.3

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 02_PRECAST_run_nnSVG.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
