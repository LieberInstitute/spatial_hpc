#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -pe local 5
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code/logs/countSpots.$TASK_ID.txt
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code/logs/countSpots.$TASK_ID.txt
#$ -m e
#$ -M madhavitippani28@gmail.com
#$ -t 1-32
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
echo "Sample id: $(cat /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code/ALLsamples.txt | awk '{print $1}' | awk "NR==${SGE_TASK_ID}") "
echo "****"

module load matlab/R2019a

toolbox='/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code'
mask=$(cat /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code/ALLsamples.txt | awk '{print $1}' | awk "NR==${SGE_TASK_ID}")
mask=$(echo ${mask} | cut -d "." -f 1)
add="_nuclei_WS_final.mat"
maskname=$mask$add

path1=$(cat /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg/code/ALLsamples.txt |  awk '{print $3}' | awk "NR==${SGE_TASK_ID}")
sample=$(echo ${mask} | cut -d "/" -f 11)

rest="/outs/spatial/scalefactors_json.json"
jsonname=$path1$sample$rest

rest1="/outs/spatial/tissue_positions_list.csv"
posname=$path1$sample$rest1

matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei('$maskname','$jsonname','$posname')"

echo "**** Job ends ****"
date
