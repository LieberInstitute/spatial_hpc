#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=80G
#$ -pe local 8
#$ -N pseudobulk_sestan
#$ -o logs/pseudobulk_sestan.txt
#$ -e logs/pseudobulk_sestan.txt
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
module load conda_R/4.3

## List current modules for reproducibility
module list

## Edit with your job command
Rscript pseudobulk_sestan.R

echo "**** Job ends ****"
date