#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_stack=256M,h_fsize=100G
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg_latest/code/logs/$TASK_ID_refineVNS.txt 
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg_latest/code/logs/$TASK_ID_refineVNS.txt
#$ -m e
#$ -M heenadivecha@gmail.com
#$ -t 1-4
#$ -tc 4

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"s
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: $(cat /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg_latest/code/refineVNS_list.txt | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

## load MATLAB
module load matlab/R2019a

## Load toolbox for VistoSeg
toolbox='/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg_latest/code'

## Read inputs from refineVNS_list.txt file
fname=$(awk 'BEGIN {FS="\t"} {print $1}' refineVNS_list.txt | awk "NR==${SGE_TASK_ID}")
M=$(awk 'BEGIN {FS="\t"} {print $2}' refineVNS_list.txt | awk "NR==${SGE_TASK_ID}")

## Run refineVNS function
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), refineVNS('$fname',$M)"

echo "**** Job ends ****"
date



