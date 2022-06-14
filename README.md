# Integration of single-nucleus and spatial transcriptomics reveals the molecular landscape of the human hippocampus

This is the repository for the spatial (U01) hippocampus (HPC) project. The README.md contains a description of files in the repository including code and data to analyze the HPC data. 

## JHPCE Internal links

* JHPCE location: `/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc`
* Use file structure similar to https://github.com/LieberInstitute/Visium_IF_AD

## Description of HPC data

This is a description of data files for this project. 

- The `/spatial_hpc/raw-data/FASTQ/` contains FASTQ files for all experiments
1. The `MiSeq` and `NovaSeq` folders has softlinks to the fastqs of slides `V10B01−085` and `V10B01−086`
2. The `2022-04-12_SPag033122` folder has softlinks to the fastqs of slides `V11A20−297`, `V11L05−333`, `V11L05−335`, `V11L05−336`, `V11U08−081`, `V11U08−084`.
- The `/spatial_hpc/raw-data/images/` contains tif files captured on CS2 for all slides 
- The `/spatial_hpc/raw-data/sample_info/` contains information of all brains, visium slides and their master excel sheets used in the study

# Description of analyses of HPC data 

## SpaceRanger
- Script to run space ranger with miseq and novaseq fastqs combined for samples `V10B01−085` and `V10B01−086` is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/01_spaceranger/spaceranger_NovaSeq.sh). 
- Script to run space ranger for all other samples is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/01_spaceranger/spaceranger_2022-04-12_SPag033122.sh). 
## Cell segmentation
### VistoSeg [here]()
### cellpose [here]()

## REDCap
- Script to extract HPC info only from the redcap form and extract all relevant (demographic/biological/rotation info etc) data to add to the spe object is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/REDCap/REDCap.R)

## Build spe
- Script to build initial raw SpatialExperiment (SPE) object from spaceranger output is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/02_build_spe/01_raw_spe.R). 
- Script to perform rotations, remove out of tissue spots and rearrange capture areas to form the HPC structure and then build basic SpatialExperiment (SPE) object is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/02_build_spe/02_basic_spe.R)

## shiny app
- Script to subset the basic_spe to make it memory effecient for shiny app is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/03_shiny_app_basic/subset.R). 
- Script to deploy the shiny app is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/03_shiny_app_basic/deploy.R). 
- Scripts for running the shiny app is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/03_shiny_app_basic/app.R). 

## QC checks
- Script to build violin plots and spot plots for all samples to show outlier spots that can be discarded from analysis is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/04_QC/qc.R)


