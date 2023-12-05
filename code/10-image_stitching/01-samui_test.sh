#!/bin/bash

#$ -cwd
#$ -N "samui_test"
#$ -o /dev/null
#$ -e /dev/null
#$ -l mf=80G,h_vmem=80G,h_fsize=50G

#SBATCH -p shared
#SBATCH --mem=80G
#SBATCH --job-name=samui_test
#SBATCH -o /dev/null
#SBATCH -e /dev/null

donor="Br6522"
mode="adjusted"

if [[ ! -z $SLURMD_NODENAME ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
    node_name=$SLURMD_NODENAME
    module_name=samui
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
    node_name=$HOSTNAME
    module_name=loopy
fi

log_path="../../processed-data/02_image_stitching/01-samui_test_${donor}_${mode}.log"

{
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Node name: ${node_name}"

module load ${module_name}/1.0.0-next.24
python 01-samui_test.py $donor $mode

echo "**** Job ends ****"
date
} > $log_path 2>&1
