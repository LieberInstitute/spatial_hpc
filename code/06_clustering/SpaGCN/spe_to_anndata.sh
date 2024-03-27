#!/bin/bash
#$ -cwd
#$ -o logs/spe_to_anndata.txt
#$ -e logs/spe_to_anndata.txt
#$ -l mf=40G,h_vmem=40G
#$ -N spe_to_anndata
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"


module load conda_R/devel
module list

Rscript spe_to_anndata.R

echo "**** Job ends ****"
date