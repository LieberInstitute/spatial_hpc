#!/bin/bash
#$ -cwd
#$ -l mem_free=40G,h_vmem=40G,h_fsize=100G
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code/logs/countSpots.$TASK_ID.txt
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code/logs/countSpots.$TASK_ID.txt
#$ -m e
#$ -M madhavitippani28@gmail.com
#$ -t 2-36
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
echo "Sample id: $(cat /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code/samples.txt | awk '{print $2}' | awk "NR==${SGE_TASK_ID}") "
echo "****"

module load matlab/R2019a

toolbox='/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code'
path='/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/Images/VistoSeg/Capture_areas/'
sample=$(cat /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code/samples.txt | awk '{print $2}' | awk "NR==${SGE_TASK_ID}")
add="_nuclei_WS.mat"
maskname=$path$sample$add

path1=$(cat /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code/samples.txt |  awk '{print $1}' | awk "NR==${SGE_TASK_ID}")
rest='outs/spatial/scalefactors_json.json'
jsonname=$path1$rest

rest1='outs/spatial/tissue_positions_list.csv'
posname=$path1$rest1

matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei('$maskname','$jsonname','$posname')"

echo "**** Job ends ****"
date
