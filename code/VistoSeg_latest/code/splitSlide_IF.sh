#!/bin/bash
#$ -l mem_free=120G,h_vmem=120G,h_fsize=100G
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg_latest/code/logs/splitSlide_IF.txt
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg_latest/code/logs/splitSlide_IF.txt
#$ -m be
#$ -M madhavitippani28@gmail.com

echo "**** Job starts ****"
date
 
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"

module load matlab/R2019a

toolbox='/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg_latest/code/'

fname1='/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/raw-data/Images/VSPG/V12D07-332_rerun/V12D07-332.mat'
fname2='/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/raw-data/Images/VSPG/V12D07-335_rerun/V12D07-335.mat'

echo "splitting V12D07-332"
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide_IF('$fname1')"
echo "splitting V12D07-335"
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide_IF('$fname2')"

echo "**** Job ends ****"
date

