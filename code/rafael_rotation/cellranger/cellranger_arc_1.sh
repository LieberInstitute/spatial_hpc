#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=8G,h_vmem=8G,h_fsize=15G
#$ -pe local 8
#$ -N round1
#$ -o logs/round1.$TASK_ID.txt
#$ -e logs/round1.$TASK_ID.txt
#$ -m e
#$ -t 1-3
#$ -tc 3

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load CellRanger
module load cellranger_arc/2.0.2

## List current modules for reproducibility
module list

## Locate file
SAMPLE=42_1
echo "Processing sample ${SAMPLE}"
date

## Run CellRanger
cellranger-arc count --id=${SAMPLE} \
    --reference=/fastscratch/myscratch/${USER}/refdata/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz \
    --libraries=/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/rafael_rotation/cellranger/libraries_1.csv \
    --localcores=8 \
    --localmem=64 \

## Move output
echo "Moving data to new location"
date
mkdir -p /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger/
mv ${SAMPLE} /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger/

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
