#!/bin/bash
#$ -cwd
#$ -l mem_free=14G,h_vmem=14G,h_fsize=112G
#$ -pe local 8
#$ -N spatialPreprocess_harmony_id
#$ -o logs/spatialPreprocess_harmony_allSamples_id.txt
#$ -e logs/spatialPreprocess_harmony_allSamples_id.txt
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
Rscript OSCApreprocess_harmony_captureArea_allSamples.R