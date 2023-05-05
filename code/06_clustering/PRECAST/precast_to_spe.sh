#!/bin/bash
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G,h_fsize=50G
#$ -N add_precast2spe
#$ -o logs/add_precast2spe.txt
#$ -e logs/add_precast2spe.txt
#$ -m e
RDATA_FILE="/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/06_clustering/PRECAST/allSamples_PRECASTObj_30.Rdata"
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
Rscript add_PRECAST_to_spe.R