#!/bin/bash -l

#$ -cwd
#$ -N "clean_sample_info"
#$ -o ../../processed-data/02_image_stitching/02-clean_sample_info.log
#$ -e ../../processed-data/02_image_stitching/02-clean_sample_info.log
#$ -l mf=5G,h_vmem=5G

#SBATCH -p shared
#SBATCH --mem=5G
#SBATCH --job-name=clean_sample_info
#SBATCH -o ../../processed-data/02_image_stitching/02-clean_sample_info.log
#SBATCH -e ../../processed-data/02_image_stitching/02-clean_sample_info.log

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

module load samui/1.0.0-next.24
python 02-clean_sample_info.py

echo "**** Job ends ****"
date
