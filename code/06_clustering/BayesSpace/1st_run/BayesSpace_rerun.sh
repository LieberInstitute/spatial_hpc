#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=200G,h_fsize=200G
#$ -N BayesSpace_rerun
#$ -o logs/BayesSpace_rerun_k17_mbkmeans.txt
#$ -e logs/BayesSpace_rerun_k17_mbkmeans.txt
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
module load conda_R/devel

## List current modules for reproducibility
module list

## Edit with your job command
Rscript BayesSpace_rerun.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
