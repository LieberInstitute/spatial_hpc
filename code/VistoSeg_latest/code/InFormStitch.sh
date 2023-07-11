#!/bin/bash
#$ -pe local 5
#$ -l mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -o /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/raw-data/Images/VSPG/logs/InFormStitch_output.txt
#$ -e /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/raw-data/Images/VSPG/logs/InFormStitch_error.txt
#$ -m be
#$ -M uma.kaipa.94@gmail.com

echo "**** Job starts ****"
date
 
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"

module load matlab/R2019a

toolbox = '/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VistoSeg_latest/code/'

path1 = '/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/raw-data/Images/VSPG/V12D07-332_0%_overlap/20230517_VSPG_HPC_Round1_Scan2_*_component_data.tif'

path2 = '/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/raw-data/Images/VSPG/V12D07-335_0%_overlap/20230518_VSPG_HPC_Round2_Scan2_*_component_data.tif'

fname1 = '/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/raw-data/Images/VSPG/V12D07-332_0%_overlap/20230517_VSPG_HPC_Round1_Scan2_[5711,29365]_component_data.tif'

fname2 = '/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/raw-data/Images/VSPG/V12D07-335_0%_overlap/20230518_VSPG_HPC_Round2_Scan2_[4867,28412]_component_data.tif'

filename1 = 'V12D07-332'

filename2 = 'V12D07-335'

matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), O1{1} = DAPI; O1{2} = 'Alexa_488'; O1{3} = 'Alexa_594'; O1{4} = 'Alexa_647'; O1{5} = 'Autofluorescence'; InFormStitch('$path1',O1,6,'$filename1')"

matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox'));O2 = extractMD('$fname2'), O2{1} = DAPI; O2{2} = 'Alexa_488'; O2{3} = 'Alexa_594'; O2{4} = 'Alexa_647'; O2{5} = 'Autofluorescence'; InFormStitch('$path2',O2,6,'$filename2')"

echo "**** Job ends ****"
date

