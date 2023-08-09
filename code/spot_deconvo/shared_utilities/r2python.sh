#!/bin/bash
#SBATCH --mem=100G
#SBATCH --job-name=r2python
#SBATCH -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/shared_utilities/logs/r2python_broad.log
#SBATCH -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/shared_utilities/logs/r2python_broad.log

#$ -l mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -N "r2python"
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/shared_utilities/logs/r2python_broad.log
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/shared_utilities/logs/r2python_broad.log

USE_SLURM=2

if [[ $USE_SLURM -eq 1 ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
fi

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Hostname: ${HOSTNAME}" 

module load conda_R/devel
Rscript /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/shared_utilities/r2python.R

echo "**** Job ends ****"
date