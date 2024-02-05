#!/bin/bash

#$ -cwd
#$ -N "samui_test"
#$ -o /dev/null
#$ -e /dev/null
#$ -l mf=80G,h_vmem=80G,h_fsize=50G

#SBATCH --mem=80G
#SBATCH --job-name=create_samui
#SBATCH -o /dev/null
#SBATCH -e /dev/null

donor="Br2720"
mode="initial"

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

log_path="logs/01-samui_test_${donor}_${mode}.log"

{
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Node name: ${node_name}"

module load samui/1.0.0-next.45
python 02-create_samui.py $donor $mode

echo "**** Job ends ****"
date
} > $log_path 2>&1
