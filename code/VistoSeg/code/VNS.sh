#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -pe local 7
#$ -o logs/VNS.$TASK_ID.txt
#$ -e logs/VNS.$TASK_ID.txt
#$ -m e
#$ -M madhavitippani28@gmail.com
#$ -t 2-40
#$ -tc 10


echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: $(cat ALLSAMPLES.txt | awk "NR==${SGE_TASK_ID}") "
echo "****"

module load matlab/R2019a

toolbox='/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code'
fname=$(cat ALLSAMPLES.txt | awk "NR==${SGE_TASK_ID}")

matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), VNS('$fname',5)"
echo "**** Job ends ****"
date
