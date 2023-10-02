#!/bin/bash
#$ -cwd
#$ -N "runRCTD"
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/RCTD/logs/03_runRCTD_IF_layer_newClass_RCTDmarkers$TASK_ID.log
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/RCTD/logs/03_runRCTD_IF_layer_newClass_RCTDmarkers$TASK_ID.log
#$ -l caracol,mf=100G,h_vmem=100G,h_fsize=100G
#$ -t 1:8
#$ -tc 4

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
Rscript 03_runRCTD.R

echo "**** Job ends ****"
date