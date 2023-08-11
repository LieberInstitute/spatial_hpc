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
sce_in <- "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/sce/sce_final.rda"
spe_in <- here("processed-data","02_build_spe","spe_nmf_final.rda")
#spg_in <- here("processed-data", "02_build_spe", "spe_nmf_final.rda")

out <- here("processed-data", "spot_deconvo", "shared_utilities")
#  Make sure output directories exist
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

###############################################################################
#   Main
###############################################################################

#   Load objects
load(sce_in, verbose = TRUE)
# spg <- readRDS(spe_IF_in)
load(spe_in, verbose = TRUE)

#-------------------------------------------------------------------------------
#   Convert snRNA-seq and spatial R objects to AnnData python objects
#-------------------------------------------------------------------------------

#   zellkonverter doesn't know how to convert the 'spatialCoords' slot. We'd
#   ultimately like the spatialCoords in the .obsm['spatial'] slot of the
#   resulting AnnDatas, which corresponds to reducedDims(spe)$spatial in R

reducedDims(spe)$spatial <- spatialCoords(spe)
#reducedDims(spg)$spatial <- spatialCoords(spg)

#   Use Ensembl gene IDs for rownames (not gene symbol)
rownames(sce) <- rowData(sce)$gene_id

#  convert all objects to Anndatas
source(here("code", "spot_deconvo", "shared_utilities","write_anndata.R"))

print("Converting objects to AnnDatas...")
write_anndata(sce, paste0(out,"/sce.h5ad"))
colData(spe)$dateImg = as.character(colData(spe)$dateImg)
write_anndata(spe, paste0(out,"/spe.h5ad"))
#write_anndata(spg, paste0(out,"/spg.h5ad"))

session_info()
