#!/bin/bash
#SBATCH --job-name=nnSVG_HPC
#SBATCH --output=logs/nnSVG_hpc.txt
#SBATCH --error=logs/nnSVG_hpc.txt
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-type=END
#SBATCH --mail-user=erik.nelson116@gmail.com # Please replace with the appropriate email address
echo "**** Job starts ****"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
## List current modules for reproducibility
module list
module load conda_R
Rscript nnSVG.R
