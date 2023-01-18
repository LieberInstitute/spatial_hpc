#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=75G,h_vmem=75G,h_fsize=100G
#$ -N feature_selection_dimred
#$ -o logs/feature_selection_dimred.txt
#$ -e logs/feature_selection_dimred.txt
#$ -m e
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.2.x
## List current modules for reproducibility
module list
## Edit with your job command
Rscript feature_selection_dimred.R 
echo "**** Job ends ****"
date