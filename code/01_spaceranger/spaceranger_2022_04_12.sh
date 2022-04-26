#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 4
#$ -N spatialHPC_spaceranger_2022-04-12
#$ -o logs/spaceranger_.$TASK_ID.txt
#$ -e logs/spaceranger_2022-04-12.$TASK_ID.txt
#$ -m e
#$ -t 1-8
#$ -tc 8

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load SpaceRanger
module load spaceranger/1.3.0

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(awk "NR==${SGE_TASK_ID}" /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/REDCap/samples.txt)
echo "Processing sample ${SAMPLE}"
date

## Get slide and area
SLIDE=$(echo ${SAMPLE} | cut -d "_" -f 1)
CAPTUREAREA=$(echo ${SAMPLE} | cut -d "_" -f 2)
echo "Slide: ${SLIDE}, capture area: ${CAPTUREAREA}"

## Find FASTQ file path
FASTQPATH2022-04-12=$(ls -d ../../raw-data/FASTQ/2022-04-12/${SAMPLE}/)
##FASTQPATHNOVASEQ=$(ls -d ../../raw-data/FASTQ/2022_04_12/${SAMPLE}/)

## Run SpaceRanger
spaceranger count \
    --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=${FASTQPATH2022-04-12}\
    --image=../../processed-data/Images/VistoSeg/Capture_areas/${SAMPLE}.tif \
    --slide=${SLIDE} \
    --area=${CAPTUREAREA} \
    --loupe-alignment=../../processed-data/Images/loupe-alignment/${SAMPLE}.json \
    --jobmode=local \
    --localcores=4 \
    --localmem=40

## Move output
echo "Moving results to new location"
date
mkdir -p ../../processed-data/spaceranger_2022-04-12/
mv ${SAMPLE} ../../processed-data/spaceranger_2022-04-12/

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
