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


## Pilot data 

- Scripts to run space ranger with miseq and novaseq samples are [here](). 
- Scripts to build initial SpatialExperiment (SPE) object are [here](). 
- Scripts for preprocessing and quality control are [here](). 
- Scripts for batch correction with harmony are [here](). 
- Scripts for unsupervised clustering with Bayespaces with `k=7` are [here](). 

## New data 

- Scripts to run space ranger with miseq and novaseq samples are [here](). 
- Scripts to extract demographic/biological data from REDCap forms instead of manually entering are [here]().  
- Scripts to reorder capture areas according to the rotations are [here](). 
- Scripts to combine all data into one final SPE object are [here](). 

## Shiny app

- Build shiny app
- Figure out rotation and flipping issues (update Vistoseg accordingly)
- Deploy shiny app
- Add here. 



