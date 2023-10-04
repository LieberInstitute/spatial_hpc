#!/bin/bash
#SBATCH --mem=50G
#SBATCH --job-name=RCTDprepare_data_HE
#SBATCH -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/RCTD/logs/01_prepare_data_HE_broad.log
#SBATCH -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/RCTD/logs/01_prepare_data_HE_broad.log
#SBATCH --array=1-36%8

#$ -cwd
#$ -N "RCTDprepare_data_IF"
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/RCTD/logs/02_prepare_myRCTD_IF_layer_deconvoMarkers$TASK_ID.log
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/RCTD/logs/02_prepare_myRCTD_IF_layer_deconvoMarkers$TASK_ID.log
#$ -l caracol,mf=50G,h_vmem=50G,h_fsize=50G
#$ -t 1-8
#$ -tc 8

USE_SLURM=2

if [[ $USE_SLURM -eq 1 ]]; then
    job_id=$SLURM_JOB_ID
    job_name=$SLURM_JOB_NAME
	task_id=$SLURM_ARRAY_TASK_ID
	hostname=$SLURMD_NODENAME
else
    job_id=$JOB_ID
    job_name=$JOB_NAME
	task_id=$SGE_TASK_ID
	hostname=$HOSTNAME
fi


echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${job_id}"
echo "Job name: ${job_name}"
echo "Hostname: ${hostname}"
echo "Task id: ${task_id}"

module load conda_R/devel
Rscript 02_prepare_myRCTD.R

echo "**** Job ends ****"
date