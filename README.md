Integration of single-nucleus and spatial transcriptomics reveals the
molecular landscape of the human hippocampus
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

This is the repository for the spatial (U01) hippocampus (HPC) project.
The README.md contains a description of files in the repository
including code and data to analyze the HPC data.

## Study design

<img src="https://research.libd.org/spatial_hpc/img/Copy%20of%20HPC%20figure%201.png" width="1000px" align="left" />

Experimental design to generate paired single-nucleus RNA-sequencing
(snRNA-seq) and spatially-resolved transcriptomics (SRT) data in the
human hippocampus. (A) Postmortem human tissue blocks from the anterior
hippocampus were dissected from 10 adult neurotypical brain donors.
Tissue blocks were scored and cryosectioned for snRNA-seq assays (red),
and placement on Visium slides (Visium H&E, black; Visium Spatial
Proteogenomics (SPG), yellow). (B) 10μm tissue sections from all ten
donors were placed onto 2-5 capture areas to include the extent of the
HPC(n=36 total capture areas), for measurement with the 10x Genomics
Visium H&E platform. (C) 10μm tissue sections from two donors were
placed onto 4 capture areas (n=8 total capture areas) for measurement
with the 10x Genomics Visium-SPG platform. (D) Tissue sections (2-4
100μm cryosections per assay) from all ten donors were collected from
the same tissue blocks for measurement with the 10x Genomics 3’ gene
expression platform . For each donor, we sorted on both and PI+NeuN+
(n=26 total snRNA-seq libraries).

## Interactive Websites

All of these interactive websites are powered by open source software,
namely:

- 🔍 [`samui`](http://dx.doi.org/10.1017/S2633903X2300017X)
- 👀 [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1)

We provide the following interactive websites, organized by dataset with
software labeled by emojis:

- Visium (n = xx)
  - 👀 TODO
  - 🔍
    [TODO](https://samuibrowser.com/from?url=data.libd.org/samuibrowser/&s=Br3942&s=Br8325&s=Br2720&s=Br2743&s=Br3942-VSPG&s=Br6423&s=Br6432&s=Br6471&s=Br6522&s=Br8325-VSPG&s=Br8492&s=Br8667)
- snRNA-seq (n = xx)
  - 👀 TODO

## JHPCE Internal links

- JHPCE location:
  `/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc`
- Use file structure similar to
  <https://github.com/LieberInstitute/Visium_IF_AD>

## Description of HPC data

This is a description of data files for this project.

- The `/spatial_hpc/raw-data/FASTQ/` contains FASTQ files for all
  experiments

1.  The `MiSeq` and `NovaSeq` folders has softlinks to the fastqs of
    slides `V10B01−085` and `V10B01−086`
2.  The `2022-04-12_SPag033122` folder has softlinks to the fastqs of
    slides `V11A20−297`, `V11L05−333`, `V11L05−335`, `V11L05−336`,
    `V11U08−081`, `V11U08−084`.

- The `/spatial_hpc/raw-data/images/` contains tif files captured on CS2
  for all slides
- The `/spatial_hpc/raw-data/sample_info/` contains information of all
  brains, visium slides and their master excel sheets used in the study

# Description of analyses of HPC data

## REDCap

- Script to extract HPC info only from the redcap form and extract all
  relevant (demographic/biological/rotation info etc) data to add to the
  spe object is
  [here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/REDCap/REDCap.R)

## 01_spaceranger

Code to align reads using 10x SpaceRanger - Script to run space ranger
with miseq and novaseq fastqs combined for samples `V10B01−085` and
`V10B01−086` is
[here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/01_spaceranger/spaceranger_NovaSeq.sh). -
Script to run space ranger for all other samples is
[here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/01_spaceranger/spaceranger_2022-04-12_SPag033122.sh).

## 02_build_spe

Code to build initial SPE, rotate/rearrange capture areas to better
reflect anatomy, and drop spots - Script to build initial raw
SpatialExperiment (SPE) object from spaceranger output is
[01_raw_spe_allSamples.R](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/02_build_spe/01_raw_spe_allSamples.R). -
Script to perform rotations and rearrange capture areas is
[02_transform_spe_allSamples.R](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/02_build_spe/02_transform_spe_allSamples.R) -
Script to remove drop spots not in tissue section and spots with zero
counts is
[03_dropSpots.R](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/02_build_spe/03_dropSpots.R).

## 03_shiny app_basic

Code to build temporary R Shiny app for initial manual annotations -
Script to subset the basic_spe to make it memory effecient for shiny app
is
[here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/03_shiny_app_basic/subset.R). -
Script to deploy the shiny app is
[here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/03_shiny_app_basic/deploy.R). -
Scripts for running the shiny app is
[here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/03_shiny_app_basic/app.R).

## 04_QC

Code to perform QC on SRT data - Script to run QC is
[here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/04_QC/qc_metrics_allSamples.R) -
Script to plot results from QC is
[here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/04_QC/qc_metrics_spotPlots_allSamples.R)

## 05_preprocess_batchCorrection

Code to preprocess data for clustering is
[here](https://github.com/LieberInstitute/spatial_hpc/blob/main/code/05_preprocess_batchCorrection/OSCApreprocess_allSamples_HE_VSPG.R)

## 06_clustering

Code to cluster SRT data. Primary method for clustering was PRECAST.
Code for formatting data, running PRECAST, and visualizing results is
[here](https://github.com/LieberInstitute/spatial_hpc/tree/main/code/06_clustering/PRECAST).

## 08_pseudobulk

Code to perform pseudobulk DE analysis is
[here](https://github.com/LieberInstitute/spatial_hpc/tree/main/code/08_pseudobulk/PRECAST)

## nnSVG

Code to identify spatially variable genes (SVGs) is
[here](https://github.com/LieberInstitute/spatial_hpc/tree/main/code/nnSVG_)

## spot_deconvo

Code for running spot-level deconvolution is
[here](https://github.com/LieberInstitute/spatial_hpc/tree/main/code/spot_deconvo)

## NMF

Code to project NMF patterns learned in paired snRNA-seq data to Visium
data and visualize results is
[here](https://github.com/LieberInstitute/spatial_hpc/tree/main/code/NMF)

## enrichment_analysis

Code for running LDSC across spatial domains, snRNA-seq cell classes,
and NMF patterns is
[here](https://github.com/LieberInstitute/spatial_hpc/tree/main/code/enrichment_analysis)

## Cell segmentation

### VistoSeg [here]()

### cellpose [here]()
