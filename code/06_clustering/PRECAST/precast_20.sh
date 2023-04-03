#!/bin/bash
#$ -cwd
#$ -l mem_free=16G,h_vmem=16G,h_fsize=16G
#$ -N plots_k_20_precast
#$ -o logs/plots_k_20_precast.txt
#$ -e logs/plots_k_20_precast.txt
#$ -m e
RDATA_FILE="/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/06_clustering/PRECAST/stitch_PRECASTObj_20.Rdata"
echo "**** HOLDING FOR RDATA OBJECT ****"
# Wait for the Rdata file to be available
while [ ! -f $RDATA_FILE ]; do
sleep 600 # Wait for 600 seconds (10 minutes) before checking again
done
echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
## List current modules for reproducibility
module list
module load conda_R
Rscript 02_PRECAST_20.R