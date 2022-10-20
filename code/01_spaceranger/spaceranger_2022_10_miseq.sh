#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 4
#$ -N spatialHPC_spaceranger_miseq_10v
#$ -o logs/spaceranger_2022_10_miseq_10v.txt
#$ -e logs/spaceranger_2022_10_miseq_10v.txt
#$ -m e
#$ -t 5
#$ -tc 5

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
SAMPLE=$(awk "NR==${SGE_TASK_ID}" samples_2022_10.txt)
echo "Processing sample ${SAMPLE}"
date

## Get slide and area
SLIDE=$(echo ${SAMPLE} | cut -d "_" -f 1)
CAPTUREAREA=$(echo ${SAMPLE} | cut -d "_" -f 2)
echo "Slide: ${SLIDE}, capture area: ${CAPTUREAREA}"

## Find FASTQ file path
FASTQPATHMISEQ=$(ls -d ../../raw-data/FASTQ/2022-10-12_MiSeq_SCP_visium/V12F14-051_A1/10v/)

## Run SpaceRanger
spaceranger count \
    --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=${FASTQPATHNOVASEQ} \
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
mkdir -p ../../processed-data/01_spaceranger/spaceranger_2022_10/V12F14-051_A1/10v/
mv ${SAMPLE} ../../processed-data/01_spaceranger/spaceranger_2022_10/V12F14-051_A1/10v/

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
