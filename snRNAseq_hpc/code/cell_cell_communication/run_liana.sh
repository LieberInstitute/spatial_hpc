#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=70G,h_vmem=70G,h_fsize=70G
#$ -N run_liana_prelim
#$ -o logs/run_liana.txt
#$ -e logs/run_liana.txt
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
module load conda_R
## List current modules for reproducibility
module list
## Edit with your job command
Rscript /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/code/cell_cell_communication/run_liana.R
echo "**** Job ends ****"
date
## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/