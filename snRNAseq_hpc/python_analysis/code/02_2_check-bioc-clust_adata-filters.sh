#!/bin/bash
#SBATCH --mem=80G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jthom338@jh.edu
#SBATCH --job-name=check_clust
#SBATCH --output=./snRNAseq_hpc/python_analysis/code/%x_%j.log

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node(s): ${SLURM_NODELIST}"
echo "Node memory requested: ${SLURM_MEM_PER_NODE}"
echo "n Tasks: ${SLURM_NTASKS}"
date
module load conda_R/devel
module list
Rscript ./snRNAseq_hpc/python_analysis/code/02_2_check-bioc-clust_adata-filters.r
