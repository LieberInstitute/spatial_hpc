#!/bin/bash
#$ -cwd
#$ -N "cart_visualize"
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/groundTruth/03_CART/logs/cart_visualize$TASK_ID.log
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/groundTruth/03_CART/logs/cart_visualize$TASK_ID.log
#$ -l caracol,mf=20G,h_vmem=20G
#$ -t 1-6
#$ -tc 6

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load cellpose/2.0
python cart_visualize.py

echo "**** Job ends ****"
date