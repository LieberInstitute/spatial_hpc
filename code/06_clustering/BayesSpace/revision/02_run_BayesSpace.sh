#!/bin/bash
#SBATCH --mem=100G
#SBATCH --job-name=run_BayesSpace
#SBATCH -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/06_clustering/BayesSpace/revision/logs/%x_%j.log


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node(s): ${SLURM_NODELIST}"
echo "Node memory requested: ${SLURM_MEM_PER_NODE}"
echo "n Tasks: ${SLURM_NTASKS}"

## Load the R module
module load conda_R/

## List current modules for reproducibility
module list

## Edit with your job command
Rscript /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/06_clustering/BayesSpace/revision/02_run_BayesSpace.r

echo "**** Job ends ****"
date

