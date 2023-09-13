setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("tidyverse"))

#  Paths
Dr <- here("processed-data","spot_deconvo","shared_utilities")

#   Load objects
# sce = readRDS(here(Dr,"sce.rds"))
sce = readRDS(here(Dr,"sce_class.rds"))
spe = readRDS(here(Dr,"spe.rds"))
#spg =  readRDS(here(Dr,"spg.rds"), verbose = TRUE)

cell_group = "broad"
cell_type_var = "broad.class"
name = "_class"

cell_group = "layer" 
cell_type_var = "cell.class"
name = "_celltype_class"

n_markers_per_type <- 25

print(paste0("Running script at ", cell_group, "-resolution."))

#-------------------------------------------------------------------------------
#   Rank marker genes
#-------------------------------------------------------------------------------

#   We won't consider genes that aren't in both the spatial objects
# keep <- (rowData(sce)$gene_id %in% rowData(spe_IF)$gene_id) &
#   (rowData(sce)$gene_id %in% rowData(spe_HE)$gene_id)

keep = (rowData(sce)$gene_id %in% rowData(spe)$gene_id)
sce <- sce[keep, ]

perc_keep <- 100 * (1 - length(which(keep)) / length(keep))
print(
  paste0(
    "Dropped ", round(perc_keep, 1), "% of potential marker genes ",
    "that were not present in all the spatial data"
  )
)

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

#-------------------------------------------------------------------------------
#   Filter out mitochondrial genes and re-rank 'rank_ratio' values. Add gene symbols to 'marker_stats' object
#-------------------------------------------------------------------------------
source(here("code", "spot_deconvo", "shared_utilities", "shared_function.R"))
#   Add gene symbol
marker_stats$symbol <- rowData(sce)$gene_name[match(marker_stats$gene, rownames(sce))]

#   Filter out mitochondrial genes
marker_stats <- marker_stats[!grepl("^MT-", marker_stats$symbol), ]

#   "Re-rank" rank_ratio, since there may be missing ranks now
cell_types = unique(marker_stats$cellType.target)
for (ct in cell_types) {
  old_ranks <- marker_stats |>
    filter(cellType.target == ct) |>
    pull(rank_ratio) |>
    sort()
  
  for (i in 1:length(which((marker_stats$cellType.target == ct)))) {
    index <- which(
      (marker_stats$cellType.target == ct) &
        (marker_stats$rank_ratio == old_ranks[i])
    )
    stopifnot(length(index) == 1)
    marker_stats[index, "rank_ratio"] <- i
  }
}

###############################################################################
#  Subset and write markers
###############################################################################
source(here("code","spot_deconvo","shared_utilities","plottingfunctions.R"))
print("Writing markers...")

write_markers(n_markers_per_type, here(Dr,paste0("markers_",cell_group,name,".txt")))
saveRDS(marker_stats, here(Dr,paste0("marker_stats_",cell_group,name,".rds")))
