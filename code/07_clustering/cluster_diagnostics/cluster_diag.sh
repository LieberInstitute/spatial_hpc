#!/bin/bash
#$ -cwd
#$ -l mem_free=24G,h_vmem=24G,h_fsize=24G
#$ -N plots_diag
#$ -o logs/plots_diag.txt
#$ -e logs/plots_diag.txt
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
Rscript cluster_diag.R