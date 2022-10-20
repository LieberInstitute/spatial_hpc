#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 4
#$ -N spatialHPC_spaceranger_miseq_20v
#$ -o logs/spaceranger_2022_10_miseq_20v.txt
#$ -e logs/spaceranger_2022_10_miseq_20v.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## load SpaceRanger
module load spaceranger/1.3.0

## List current modules for reproducibility
module list

## Locate file
SAMPLE=V12F14-051_A1
echo "Processing sample ${SAMPLE}"
date

## Get slide and area
SLIDE=V12F14-051
CAPTUREAREA=A1
echo "Slide: ${SLIDE}, capture area: ${CAPTUREAREA}"

## Find FASTQ file path
FASTQPATHMISEQ=$(ls -d ../../raw-data/FASTQ/2022-10-12_MiSeq_SCP_visium/V12F14-051_A1/20v/)

## Run SpaceRanger
spaceranger count \
    --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=${FASTQPATHMISEQ} \
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
mkdir -p ../../processed-data/01_spaceranger/spaceranger_2022_10/V12F14-051_A1/20v/
mv ${SAMPLE} ../../processed-data/01_spaceranger/spaceranger_2022_10/V12F14-051_A1/20v/

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
