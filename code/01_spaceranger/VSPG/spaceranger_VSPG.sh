#!/bin/bash
#$ -cwd
#$ -l mem_free=8G,h_vmem=8G,h_fsize=100G
#$ -pe local 8
#$ -N sparanger_HPC_VSPG
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/01_spaceranger/VSPG/logs/sparanger_HPC_VSPG.$TASK_ID.txt
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/01_spaceranger/VSPG/logs/sparanger_HPC_VSPG.$TASK_ID.txt
#$ -m e
#$ -M madhavitippani28@gmail.com
#$ -t 1-8
#$ -tc 4
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load SpaceRanger
module load spaceranger/2.0.0

## List current modules for reproducibility
module list
