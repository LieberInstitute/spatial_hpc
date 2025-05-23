#!/bin/bash
#$ -cwd
#$ -l mem_free=72G,h_vmem=72G,h_fsize=72G
#$ -N NMF_firstGo
#$ -o logs/NMF_firstGo.txt
#$ -e logs/NMF_firstGo.txt
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
Rscript RcppML_nmf.R