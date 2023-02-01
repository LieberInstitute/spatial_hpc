#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_stack=256M,h_fsize=100G
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg_latest/code/logs/$TASK_ID.txt
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg_latest/code/logs/$TASK_ID.txt
#$ -m e
#$ -M heenadivecha@gmail.com
#$ -t 1
#$ -tc 1

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: $(cat /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg_latest/code/VNS_list.txt | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

module load matlab/R2019a


toolbox='/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg_latest/code'
fname=$(cat /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg_latest/code/VNS_list.txt | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")

matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), VNS('$fname',5)"

echo "**** Job ends ****"
date



