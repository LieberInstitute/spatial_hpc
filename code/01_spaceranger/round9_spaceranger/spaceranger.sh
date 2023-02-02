#!/bin/bash
#$ -cwd
#$ -l mem_free=8G,h_vmem=8G,h_fsize=100G
#$ -pe local 8
#$ -N spatial_hpc
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/01_spaceranger/round9_spaceranger/logs/round9_spaceranger.$TASK_ID.txt
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/01_spaceranger/round9_spaceranger/logs/round9_spaceranger.$TASK_ID.txt
#$ -m e
#$ -M heenadivecha@gmail.com
#$ -t 1-4
#$ -tc 10

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load SpaceRanger
module load spaceranger/1.3.1

## List current modules for reproducibility
module list

## Read parameters
SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' spaceranger_parameters.txt | awk "NR==${SGE_TASK_ID}")
SLIDE=$(awk 'BEGIN {FS="\t"} {print $2}' spaceranger_parameters.txt | awk "NR==${SGE_TASK_ID}")
CAPTUREAREA=$(awk 'BEGIN {FS="\t"} {print $3}' spaceranger_parameters.txt | awk "NR==${SGE_TASK_ID}")
IMAGEPATH=$(awk 'BEGIN {FS="\t"} {print $4}' spaceranger_parameters.txt | awk "NR==${SGE_TASK_ID}")
LOUPEPATH=$(awk 'BEGIN {FS="\t"} {print $5}' spaceranger_parameters.txt | awk "NR==${SGE_TASK_ID}")
FASTQPATH=$(awk 'BEGIN {FS="\t"} {print $6}' spaceranger_parameters.txt | awk "NR==${SGE_TASK_ID}")

echo "Processing sample ${SAMPLE} from slide ${SLIDE} and capture area ${CAPTUREAREA} with image ${IMAGEPATH} and aligned with ${LOUPEPATH} with FASTQs: ${FASTQPATH}"
date

## For keeping track of dates of the input files
ls -lh ${IMAGEPATH}
ls -lh ${LOUPEPATH}

## Hank from 10x Genomics recommended setting this environment
export NUMBA_NUM_THREADS=1

## Run SpaceRanger
spaceranger count \
    --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=${FASTQPATH} \
    --image=${IMAGEPATH} \
    --slide=${SLIDE} \
    --area=${CAPTUREAREA} \
    --loupe-alignment=${LOUPEPATH} \
    --jobmode=local \
    --localcores=8 \
    --localmem=64

## Move output
echo "Moving results to new location"
date
mkdir -p ../../../processed-data/01_spaceranger/spaceranger_2023-01-31_round9
mv ${SAMPLE} ../../../processed-data/01_spaceranger/spaceranger_2023-01-31_round9

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/