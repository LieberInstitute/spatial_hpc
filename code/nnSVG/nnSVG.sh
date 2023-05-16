#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=240G
#$ -pe local 12
#$ -N nnSVG_HPC
#$ -o logs/nnSVG_hpc.txt
#$ -e logs/nnSVG_hpc.txt
#$ -m e
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
## List current modules for reproducibility
module list
module load conda_R
Rscript nnSVG.R