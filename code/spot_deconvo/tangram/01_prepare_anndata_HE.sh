#!/bin/bash
#$ -cwd
#$ -N "prepare_anndatas_HE"
#$ -o ../../../code/spot_deconvo/tangram/logs/01_prepare_anndatas_HE_layer_$TASK_ID.log
#$ -e ../../../code/spot_deconvo/tangram/logs/01_prepare_anndatas_HE_layer_$TASK_ID.log
#$ -l mf=120G,h_vmem=120G,h_fsize=50G
#$ -t 1-36
#$ -tc 4

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load tangram/1.0.2
python 02_prepare_anndats_He.py

echo "**** Job ends ****"
date