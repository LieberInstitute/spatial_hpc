# Description of analyses of HPC data

## REDCap
- Script to extract HPC info only from the redcap form and extract all relevant (demographic/biological/rotation info etc) data to add to the spe object is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/REDCap/REDCap.R)

## 01_spaceranger
Code to align reads using 10x SpaceRanger
- Script to run space ranger with miseq and novaseq fastqs combined for samples `V10B01−085` and `V10B01−086` is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/01_spaceranger/spaceranger_NovaSeq.sh).
- Script to run space ranger for all other samples is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/01_spaceranger/spaceranger_2022-04-12_SPag033122.sh).

## 02_build_spe
Code to build initial SPE, rotate/rearrange capture areas to better reflect anatomy, and drop spots
- Script to build initial raw SpatialExperiment (SPE) object from spaceranger output is [01_raw_spe_allSamples.R](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/02_build_spe/01_raw_spe_allSamples.R).
- Script to perform rotations and rearrange capture areas is [02_transform_spe_allSamples.R](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/02_build_spe/02_transform_spe_allSamples.R)
- Script to remove drop spots not in tissue section and spots with zero counts is [03_dropSpots.R](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/02_build_spe/03_dropSpots.R).

## 03_shiny app_basic
Code to build temporary R Shiny app for initial manual annotations
- Script to subset the basic_spe to make it memory effecient for shiny app is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/03_shiny_app_basic/subset.R).
- Script to deploy the shiny app is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/03_shiny_app_basic/deploy.R).
- Scripts for running the shiny app is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/03_shiny_app_basic/app.R).

## 04_QC
Code to perform QC on SRT data
- Script to run QC is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/04_QC/qc_metrics_allSamples.R)
- Script to plot results from QC is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/04_QC/qc_metrics_spotPlots_allSamples.R)

## 05_preprocess_batchCorrection
Code to preprocess data for clustering is [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/05_preprocess_batchCorrection/OSCApreprocess_allSamples_HE_VSPG.R)

## 06_clustering
Code to cluster SRT data. Primary method for clustering was PRECAST. Code for formatting data, running PRECAST, and visualizing results is [here](https://github.com/LieberInstitute/spatial_hpc/tree/main/code/06_clustering/PRECAST).

## 08_pseudobulk
Code to perform pseudobulk DE analysis is [here](https://github.com/LieberInstitute/spatial_hpc/tree/main/code/08_pseudobulk/PRECAST)

## nnSVG
Code to identify spatially variable genes (SVGs) is [here](https://github.com/LieberInstitute/spatial_hpc/tree/main/code/nnSVG_)

## spot_deconvo
Code for running spot-level deconvolution is [here](https://github.com/LieberInstitute/spatial_hpc/tree/main/code/spot_deconvo)

## NMF
Code to project NMF patterns learned in paired snRNA-seq data to Visium data and visualize results is [here](https://github.com/LieberInstitute/spatial_hpc/tree/main/code/NMF)

## enrichment_analysis
Code for running LDSC across spatial domains, snRNA-seq cell classes, and NMF patterns is [here](https://github.com/LieberInstitute/spatial_hpc/tree/main/code/enrichment_analysis)

## Cell segmentation
### VistoSeg [here]()
### cellpose [here]()

