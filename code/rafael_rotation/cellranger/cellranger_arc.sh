#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=16G,h_vmem=16G,h_fsize=100G
#$ -pe local 4
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
SAMPLE=$(awk "NR==${SGE_TASK_ID}" ${JOB_NAME}.txt)
echo "Processing sample ${SAMPLE}"
date

## Run CellRanger
cellranger-atac count --id
=${SAMPLE} \
    --reference=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-11-22_KMay110521_ATAC \ ## reference human genome
    --fastqs=/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/raw-data/FASTQ/${SAMPLE} \
    --sample=${SAMPLE} \ ## there is no sample argument, neither fastqs, it is used though the ids
    --jobmode=local \
    --localcores=8 \
    --localmem=64 \
    --include-introns


                        --fastqs=/home/jdoe/runs/HAWT7ADXX/outs/fastq_path \

## Move output
echo "Moving data to new location"
date
mkdir -p /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger/
mv ${SAMPLE} /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger/

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
