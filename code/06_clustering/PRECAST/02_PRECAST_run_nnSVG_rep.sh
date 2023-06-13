#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=80G
#$ -pe local 8
#$ -N rep_nnSVG_PRECAST_k19-25
#$ -o logs/nnsvg_rep_PRECAST_k.$TASK_ID.txt
#$ -e logs/nnsvg_rep_PRECAST_k.$TASK_ID.txt
#$ -m e
#$ -t 16-20
#$ -tc 3

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
Rscript 02_PRECAST_run_nnSVG_rep.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
