#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G,bluejay
#$ -pe local 4
#$ -N round0_1
#$ -o logs/1c-scp.txt
#$ -e logs/1c-scp.txt
#$ -m e


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"
echo "Task id: 1c-scp"
echo "****"

##load cellranger 7
module load cellranger/7.0.0

## List current modules for reproducibility
module list


cellranger count --id=1c-scp \
     --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
		 --fastqs=/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/raw-data/FASTQ/1c-scp \
		 --sample=1c-scp \
		 --jobmode=local \
		 --localcores=4 \
		 --localmem=20 \
     