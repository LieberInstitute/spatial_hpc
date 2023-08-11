setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("HDF5Array"))

sce_in <- "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/sce/sce_final.rda"
spe_in <- here("processed-data","02_build_spe","spe_nmf_final.rda")
#spg_in <- here("processed-data", "02_build_spe", "spe_nmf_final.rda")

load(sce_in, verbose = TRUE)
rownames(sce) <- rowData(sce)$gene_id

# spg <- readRDS(spe_IF_in)
load(spe_in, verbose = TRUE)
reducedDims(spe)$spatial <- spatialCoords(spe)
#reducedDims(spg)$spatial <- spatialCoords(spg)

#cell_group = "broad" 
 cell_group = "layer" 
 
out = here("processed-data","spot_deconvo","shared_utilities",paste0("marker_stats_",cell_group,".rds"))

n_markers_per_type <- 25

#   Define variables related to cell_group

if (cell_group == "broad") {
  cell_type_var <- "broad.type"
} else {
  cell_type_var <- "cell.type"
}


if (cell_group == "broad") {
  cell_types <- c('ExcN', 'InhN', 'Glia', 'Immune', 'CSF', 'Vascular')
  colors_col <- "cell_type_colors_broad"
  cell_column <- "broad.type"
  cell_type_nrow <- 2
} else {
  cell_types <- c('GC', 'CA2-4', 'CA1', 'ProS/Sub', 'L2/3', 'L5', 'L6/6b', 'HATA/AHi',
                  'Thal', 'Cajal', 'GABA', 'Oligo', 'Astro', 'OPC', 'Micro/Macro/T',
                  'Ependy', 'Choroid', 'Vascular')
  colors_col <- "cell_type_colors_layer"
  cell_column <- "cell.type"
  cell_type_nrow <- 3
}

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

#   Change "/" to "_" for layer-level cell types
marker_stats <- marker_stats |>
  mutate(
    cellType.target = gsub("/", "_", cellType.target),
    cellType = gsub("/", "_", cellType)
  )
stopifnot(
  identical(sort(unique(marker_stats$cellType.target)), sort(cell_types))
)

#   "Re-rank" rank_ratio, since there may be missing ranks now
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
write_markers(n_markers_per_type, marker_out)

saveRDS(marker_stats, out)
