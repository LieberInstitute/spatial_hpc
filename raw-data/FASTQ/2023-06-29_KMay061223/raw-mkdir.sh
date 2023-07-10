#!/bin/bash
#$ -cwd
#$ -t 1-8
#$ -tc 8

SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' raw-mkdir.txt | awk "NR==${SGE_TASK_ID}")

mkdir -p ./logs/
mkdir -p /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/raw-data/FASTQ/2023-06-29_KMay061223/${SAMPLE}/

mv ./raw-mkdir.sh.e* ./logs
mv ./raw-mkdir.sh.o* ./logs