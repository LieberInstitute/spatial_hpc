#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 8
#$ -N OSCApreprocess_harmony_captureArea
#$ -o logs/OSCApreprocess_harmony_captureArea.txt
#$ -e logs/OSCApreprocess_harmony_captureArea.txt
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

module load conda_R/devel
Rscript OSCApreprocess_harmony_captureArea.R