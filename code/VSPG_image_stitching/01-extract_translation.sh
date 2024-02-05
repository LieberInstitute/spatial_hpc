#!/bin/bash -l

#$ -cwd
#$ -N "clean_sample_info"
#$ -o logs/02-clean_sample_info.log
#$ -e logs/02-clean_sample_info.log
#$ -l mf=5G,h_vmem=5G


#SBATCH --mem=5G
#SBATCH --job-name=extract_translation
#SBATCH -o logs/01-extract_translation.py.log


set -e

if [[ ! -z $SLURMD_NODENAME ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
    node_name=$SLURMD_NODENAME
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
    node_name=$HOSTNAME
fi

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Node name: ${node_name}"

module list

module load samui/1.0.0-next.45
python 01-extract_translation.py

echo "**** Job ends ****"
date
