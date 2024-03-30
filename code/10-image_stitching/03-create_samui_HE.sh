#!/bin/bash

#$ -cwd
#$ -N "samui_test"
#$ -o /dev/null
#$ -e /dev/null
#$ -l mf=80G,h_vmem=80G,h_fsize=50G

#SBATCH --mem=90G
#SBATCH --job-name=03-create_samui_HE
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-7%4

#donor="Br2743"
#enter donor number as bash argument
#mode="initial"
#donor="Br${1}"

donor=$(awk "NR==${SLURM_ARRAY_TASK_ID}" sample_list.txt)
echo "Processing sample ${donor}"
date

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

#log_path="logs/01-samui_test_${donor}.log"
log_path="logs/03-create_samui_HE_${donor}.log"

{
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Node name: ${node_name}"

module load samui/1.0.0-next.49
python 03-create_samui_HE.py $donor

echo "**** Job ends ****"
date
} > $log_path 2>&1
