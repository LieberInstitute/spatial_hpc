#!/bin/bash
#SBATCH --mem=120G
#SBATCH --job-name=prepare_anndatas_HE
#SBATCH -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/tangram/logs/01_prepare_anndatas_HE_broad_$SLURM_ARRAY_TASK_ID.log
#SBATCH -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/tangram/logs/01_prepare_anndatas_HE_broad_$SLURM_ARRAY_TASK_ID.log
#SBATCH --array=1-36%4

#$ -cwd
#$ -N "prepare_anndatas_HE"
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/tangram/logs/01_prepare_anndatas_HE_broad_class$TASK_ID.log
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/spot_deconvo/tangram/logs/01_prepare_anndatas_HE_broad_class$TASK_ID.log
#$ -l caracol,mf=120G,h_vmem=120G,h_fsize=50G
#$ -t 1-36
#$ -tc 2

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

module load tangram/1.0.2
python 01_prepare_anndata_HE.py

echo "**** Job ends ****"
date