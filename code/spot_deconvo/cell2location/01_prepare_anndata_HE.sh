#!/bin/bash
#$ -cwd
#$ -N "prepare_anndata_nonIF"
#$ -o ../../../processed-data/spot_deconvo/cell2location/01-prepare_anndata_HE_broad.log
#$ -e ../../../processed-data/spot_deconvo/cell2location/01-prepare_anndata_HE_broad.log
#$ -l caracol,mf=80G,h_vmem=80G,h_fsize=50G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load cell2location/0.8a0
python 01-prepare_anndata_HE.py

echo "**** Job ends ****"
date