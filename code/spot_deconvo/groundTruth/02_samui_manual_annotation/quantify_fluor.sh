#!/bin/bash
#$ -cwd
#$ -N "quantify_fluor"
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/groundTruth/02_samui_manual_annotation/logs/quantify_fluor_$TASK_ID.log
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/groundTruth/02_samui_manual_annotation/logs/quantify_fluor_$TASK_ID.log
#$ -l mf=20G,h_vmem=20G
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

module load cellpose/2.0
python quantify_fluor.py

echo "**** Job ends ****"
date