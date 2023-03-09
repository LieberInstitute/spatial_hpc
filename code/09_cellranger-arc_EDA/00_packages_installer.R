########################################################################
## 1st EDA with scCellRanger-ARC
## Author. CSC / LIEBER
## Date. March 3rd, 2023
## LM. March 7th, 2023
## Reference Link. 
## Use this first https://stuartlab.org/signac/articles/pbmc_multiomic.html 
## 
## PACKAGES NEED IT
## 
########################################################################

setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
if (!requireNamespace("Signac", quietly = TRUE))
  install.packages("Signac")
## R package for the analysis of single-cell chromatin data, including scATAC-seq, single-cell targeted tagmentation methods such as scCUT&Tag and scNTT-seq, ## and multimodal datasets that jointly measure chromatin state alongside other modalities.
if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages("Seurat")
if (!requireNamespace("SeuratDisk", quietly = TRUE))
  install.packages("SeuratDisk")
if (!requireNamespace("qlcMatrix", quietly = TRUE))    # For linking peaks to genes
  install.packages("qlcMatrix")
## Seurat is an R package designed for QC, analysis, and exploration of single-cell RNA-seq data. ## Seurat aims to enable users to identify and interpret sources of heterogeneity from single-cell transcriptomic measurements, and to integrate diverse types of single-cell data.
if (!requireNamespace("hdf5r", quietly = TRUE))
  install.packages("hdf5r")
# hdf5r is an R interface to the HDF5 library. Library has a file format to manage, process, and store your heterogeneous data.
# First need to be installed accordingly with the OS. E.g. $ sudo apt-get install libhdf5-dev
if (!require("BiocManager", quietly = TRUE))
  #  install.packages("BiocManager")
  BiocManager::install(version = "3.16")
if (!require("GenomeInfoDb", quietly = TRUE))
  BiocManager::install("GenomeInfoDb")
## GenomeInfoDb is a collection of utilities for manipulating chromosome names, including modifying them to follow a particular naming style.
## Exposes an annotation databases for human gene and protein annotations defined in Ensembl.
#if (!require("EnsDb.Hsapiens.v75", quietly = TRUE))
#  BiocManager::install("EnsDb.Hsapiens.v75")
if (!require("EnsDb.Hsapiens.v86", quietly = TRUE))
  BiocManager::install("EnsDb.Hsapiens.v86")
if (!require("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE))
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
# Data visualization package for statistics. # The easiest way to get ggplot2 is to install the whole tidyverse: install.packages("tidyverse")
if (!requireNamespace("patchwork", quietly = TRUE))
  install.packages('patchwork')
# Patchwork combine separate ggplots into the same graphic. 
