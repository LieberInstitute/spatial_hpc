#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 4
#$ -N BayesSpace
#$ -o /code/pilot_data_checks/BayesSpace.txt
#$ -e /code/pilot_data_checks/BayesSpace.txt
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

module load conda_R/4.1.x
Rscript /code/pilot_data_checks/03_BayesSpace.R