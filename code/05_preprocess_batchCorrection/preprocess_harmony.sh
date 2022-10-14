#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=12G,h_vmem=12G,h_fsize=100G
#$ -pe local 8
#$ -N spatialPreprocess_harmony
#$ -o logs/spatialPreprocess_harmony.txt
#$ -e logs/spatialPreprocess_harmony.txt
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
Rscript spatialPreprocess_harmony.R