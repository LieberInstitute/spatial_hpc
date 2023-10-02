#!/bin/bash
#$ -cwd
#$ -N "registrationHE"
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/cell2location/02_registration_IF_layer_class.log
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/cell2location/02_registration_IF_layer_class.log
#$ -l gpu,mf=300G,h_vmem=300G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

###############################################################################
#   Dynamically select a GPU based on availability
###############################################################################

# USAGE_CUTOFF=10
# NUM_GPUS=1

# avail_gpus=$(
#     nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader |
#     cut -d " " -f 1 | awk -v usage="$USAGE_CUTOFF" '$1 < usage {print NR - 1}'
# )

#  Simply exit with an error if there are no GPUs left
# if [[ -z $avail_gpus ]]; then
#     echo "No GPUs are available."
#     exit 1
# fi

# export CUDA_VISIBLE_DEVICES=$(
#     echo "$avail_gpus" | head -n $NUM_GPUS | paste -sd ","
# )

# echo "Chose GPU(s): $CUDA_VISIBLE_DEVICES"

###############################################################################
#   Submit the python script
###############################################################################

module load cell2location/0.8a0
python 02_registration_HE.py

echo "**** Job ends ****"
date