#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N barcode_ranks_empty_drops_qc
#$ -o logs/02_barcode_ranks_empty_drops_qc.$TASK_ID.txt
#$ -e logs/02_barcode_ranks_empty_drops_qc.$TASK_ID.txt
#$ -m e
#$ -t 1-26
#$ -tc 5

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 02_barcode_ranks_empty_drops_qc.R $SGE_TASK_ID

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
