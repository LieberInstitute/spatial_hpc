#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -N cluster
#$ -o logs/cluster.txt
#$ -e logs/cluster.txt
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
Rscript cluster.R 
echo "**** Job ends ****"
date