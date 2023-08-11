setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("tidyverse"))

#-------------------------------------------------------------------------------
#   Rank marker genes
#-------------------------------------------------------------------------------

#   We won't consider genes that aren't in both the spatial objects
# keep <- (rowData(sce)$gene_id %in% rowData(spe_IF)$gene_id) &
#   (rowData(sce)$gene_id %in% rowData(spe_HE)$gene_id)

sce_in <- "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/sce/sce_final.rda"
spe_HE_in <- here("processed-data","02_build_spe","spe_nmf_final.rda")
#spe_IF_in <- here("processed-data", "02_build_spe", "spe_nmf_final.rda")

load(sce_in, verbose = TRUE)
# spg <- readRDS(spe_IF_in)
load(spe_HE_in, verbose = TRUE)

cell_group = "broad" 
# cell_group = "layer" 

out = here("processed-data","spot_deconvo","shared_utilities",paste0("marker_stats_",cell_group,".rds"))

reducedDims(spe)$spatial <- spatialCoords(spe)
#reducedDims(spg)$spatial <- spatialCoords(spg)

#   Use Ensembl gene IDs for rownames (not gene symbol)
rownames(sce) <- rowData(sce)$gene_id

keep = (rowData(sce)$gene_id %in% rowData(spe)$gene_id)
sce <- sce[keep, ]

perc_keep <- 100 * (1 - length(which(keep)) / length(keep))
print(
  paste0(
    "Dropped ", round(perc_keep, 1), "% of potential marker genes ",
    "that were not present in all the spatial data"
  )
)

if (cell_group == "broad") {
  cell_type_var <- "broad.type"
} else {
  cell_type_var <- "cell.type"
}


print("Running getMeanRatio2 and findMarkers_1vAll to rank genes as markers...")
marker_stats <- get_mean_ratio2(
  sce,
  cellType_col = cell_type_var, assay_name = "logcounts"
)
marker_stats_1vall <- findMarkers_1vAll(
  sce,
  cellType_col = cell_type_var, assay_name = "logcounts",
  mod = "~brnum"
)
marker_stats <- left_join(
  marker_stats, marker_stats_1vall,
  by = c("gene", "cellType.target")
)

saveRDS(marker_stats, out)
