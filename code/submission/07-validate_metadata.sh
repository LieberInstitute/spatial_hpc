#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=3G
#SBATCH --job-name=07-validate_metadata
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -o ../../processed-data/nda-submission/logs/07-validate_metadata.txt
#SBATCH -e ../../processed-data/nda-submission/logs/07-validate_metadata.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load nda-tools/0.3.0
module list

repo_dir=$(git rev-parse --show-toplevel)
meta_dir=$repo_dir/processed-data/nda-submission

#   The tool doesn't appear to have a (documented) way to redirect results to
#   a shared location (--log_dir doesn't specify the main results output
#   folder). I interactively examined results after running this script.

vtcmd $meta_dir/rna_seq.csv
vtcmd $meta_dir/visium_image.csv
vtcmd $meta_dir/genomics_subject.csv

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/
