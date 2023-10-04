#!/bin/bash
#$ -cwd
#$ -N "cart"
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/groundTruth/03_CART/logs/cart_$TASK_ID.log
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/groundTruth/03_CART/logs/cart_$TASK_ID.log
#$ -l caracol,mf=10G,h_vmem=10G
#$ -t 7-8
#$ -tc 2

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load cellpose/2.0
python cart_test.py

echo "**** Job ends ****"
date