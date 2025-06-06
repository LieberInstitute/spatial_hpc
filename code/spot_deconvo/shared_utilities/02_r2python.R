#   Convert IF SpatialExperiment, non-IF SpatialExperiment, and snRNA-seq
#   SingleCellExperiment to AnnDatas. Save a copy of the modified SCE as an R
#   object as well. This is a processing step shared by, and prior to, the
#   different deconvolution softwares (tangram, cell2location, SPOTlight).
#   Finally, run getMeanRatio2 to rank genes as markers, and save the resulting
#   object.
setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages(library("basilisk"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("zellkonverter"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("tidyverse"))

#  Paths
Dr <- here("processed-data","spot_deconvo","shared_utilities")

###############################################################################
#   Main
###############################################################################

#   Load objects
#sce = readRDS(here(Dr,"sce.rds"), verbose = TRUE)

sce = readRDS(here(Dr,"sce_class1_noHATAGABAAmy.rds"))
spe = readRDS(here(Dr,"spe.rds"))
spg = readRDS(here(Dr,"spg.rds"))

#-------------------------------------------------------------------------------
#   Convert snRNA-seq and spatial R objects to AnnData python objects
#-------------------------------------------------------------------------------

#  convert all objects to Anndatas
source(here("code", "spot_deconvo", "shared_utilities","write_anndata.R"))

print("Converting objects to AnnDatas...")
#write_anndata(sce, here(Dr,"sce.h5ad"))
write_anndata(sce, here(Dr,"sce_class1_noHATAGABAAmy.h5ad"))
colData(spe)$dateImg = as.character(colData(spe)$dateImg)
write_anndata(spe, here(Dr,"spe.h5ad"))
colData(spg)$dateImg = as.character(colData(spg)$dateImg)
write_anndata(spg, paste0(Dr,"spg.h5ad"))

session_info()
