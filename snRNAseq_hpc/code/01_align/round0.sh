#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G,bluejay
#$ -pe local 4
#$ -N round0
#$ -o logs/round0.$TASK_ID.txt
#$ -e logs/round0.$TASK_ID.txt
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
echo "****"
echo "Task id: ${SGE_TASK_ID}"
echo "****"

##load cellranger 7
module load cellranger/7.0.0

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(awk "NR==${SGE_TASK_ID}" ${JOB_NAME}.txt)
echo "Processing sample ${SAMPLE}"
date


cellranger count --id=${SAMPLE} \
     --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
		 --fastqs=/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/raw-data/FASTQ/${SAMPLE} \
		 --sample=${SAMPLE} \
		 --jobmode=local \
		 --localcores=4 \
		 --localmem=40 \
     