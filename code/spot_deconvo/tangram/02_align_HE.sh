#!/bin/bash
#$ -cwd
#$ -N "align_tangram_HE"
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/tangram/logs/02_align_tangram_HE_layer_class$TASK_ID.log
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/tangram/logs/02_align_tangram_HE_layer_class$TASK_ID.log
#$ -l caracol,mf=210G,h_vmem=210G
#$ -t 1-8
#$ -tc 1

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

###############################################################################
#   Dynamically select a GPU based on availability
###############################################################################

#USAGE_CUTOFF=10
#NUM_GPUS=1

#avail_gpus=$(
#    nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader |
#    cut -d " " -f 1 | awk -v usage="$USAGE_CUTOFF" '$1 < usage {print NR - 1}'
#)

#  Simply exit with an error if there are no GPUs left
#if [[ -z $avail_gpus ]]; then
#    echo "No GPUs are available."
#    exit 1
#fi

#export CUDA_VISIBLE_DEVICES=$(
#    echo "$avail_gpus" | head -n $NUM_GPUS | paste -sd ","
#)

#echo "Chose GPU(s): $CUDA_VISIBLE_DEVICES"

###############################################################################
#   Submit the python script
###############################################################################

module load tangram/1.0.2

python 02_align_HE.py

echo "**** Job ends ****"
date
