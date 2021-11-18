#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 8
#$ -N spaceranger
#$ -o logs/spaceranger.$TASK_ID.txt
#$ -e logs/spaceranger.$TASK_ID.txt
#$ -m e
#$ -t 1-10
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
SAMPLE=$(awk "NR==${SGE_TASK_ID}" ../01_bamtofastq/samples.txt)
echo "Processing sample ${SAMPLE}"
date

## Get slide and area
SLIDEPART1=$(echo ${SAMPLE} | cut -c1-6)
SLIDEPART2=$(echo ${SAMPLE} | cut -c7-9)
SLIDE="${SLIDEPART1}-${SLIDEPART2}"
CAPTUREAREA=$(echo ${SAMPLE} | cut -c11-12)

## Get VIF part
if [ ${SLIDE} == "V10A27-004" ]
then
    SLIDEVIF="VIFAD1"
elif [ ${SLIDE} == "V10A27-106" ]
then
    SLIDEVIF="VIFAD2"
elif [ ${SLIDE} == "V10T31-036" ]
then
    SLIDEVIF="VIFAD3"
else
    echo "Unsupported slide ${SLIDE}."
    exit 1
fi
echo "Slide: ${SLIDE}, capture area: ${CAPTUREAREA}, VIF: ${SLIDEVIF}"

## Find FASTQ file path
FASTQPATH=$(ls -d ../../raw-data/FASTQ/spaceranger_our_alignments_nocr11/${SAMPLE}/*/)

## Run SpaceRanger
spaceranger count \
    --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=${FASTQPATH} \
    --darkimage=../../processed-data/Images/loupe_alignment/${SLIDEVIF}_${SLIDE}_${CAPTUREAREA}.tif \
    --slide=${SLIDE} \
    --area=${CAPTUREAREA} \
    --loupe-alignment=../../processed-data/Images/loupe_alignment/${SLIDEVIF}_${SLIDE}-${CAPTUREAREA}.json \
    --jobmode=local \
    --localcores=8 \
    --localmem=80

## Move output
echo "Moving results to new location"
date
mkdir -p ../../processed-data/spaceranger/
mv ${SAMPLE} ../../processed-data/spaceranger/

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
