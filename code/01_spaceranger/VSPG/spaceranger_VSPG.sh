#!/bin/bash
#$ -cwd
#$ -l mem_free=8G,h_vmem=8G,h_fsize=100G
#$ -pe local 8
#$ -N sparanger_HPC_VSPG
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/01_spaceranger/VSPG/logs/sparanger_HPC_VSPG.$TASK_ID.txt
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/01_spaceranger/VSPG/logs/sparanger_HPC_VSPG.$TASK_ID.txt
#$ -m e
#$ -M madhavitippani28@gmail.com
#$ -t 1-8
#$ -tc 4
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load SpaceRanger
module load spaceranger/2.0.0

## List current modules for reproducibility
module list

## Read parameters
BRAIN=$(awk 'BEGIN {FS="\t"} {print $1}' parameters.txt | awk "NR==${SGE_TASK_ID}")
SLIDE=$(awk 'BEGIN {FS="\t"} {print $2}' parameters.txt | awk "NR==${SGE_TASK_ID}")
CAPTUREAREA=$(awk 'BEGIN {FS="\t"} {print $3}' parameters.txt | awk "NR==${SGE_TASK_ID}")
IMAGEPATH=$(awk 'BEGIN {FS="\t"} {print $4}' parameters.txt | awk "NR==${SGE_TASK_ID}")
LOUPEPATH=$(awk 'BEGIN {FS="\t"} {print $5}' parameters.txt | awk "NR==${SGE_TASK_ID}")
FASTQPATH=$(awk 'BEGIN {FS="\t"} {print $6}' parameters.txt | awk "NR==${SGE_TASK_ID}")
SAMPLE=$(paste <(echo ${SLIDE}) <(echo "_") <(echo ${CAPTUREAREA}) <(echo "_") <(echo ${BRAIN})  -d '')

echo "Processing sample ${BRAIN} from slide ${SLIDE} and capture area ${CAPTUREAREA} with image ${IMAGEPATH} and aligned with ${LOUPEPATH} with FASTQs: ${FASTQPATH}"
date

## For keeping track of dates of the input files
ls -lh ${IMAGEPATH}
ls -lh ${LOUPEPATH}