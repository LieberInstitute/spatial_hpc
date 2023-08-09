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

cell_group <- "broad" # "broad" or "layer"

#  Paths
sce_in <- "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/sce/sce_final.rda"
#spe_IF_in <- here("processed-data", "02_build_spe", "spe_nmf_final.rda")
spe_HE_in <- here("processed-data","02_build_spe","spe_nmf_final.rda")

sce_out <- here("processed-data", "spot_deconvo", "shared_utilities", paste0("sce_", cell_group, ".h5ad"))
sce_r_out <- here("processed-data", "spot_deconvo", "shared_utilities", paste0("sce_", cell_group, ".rds"))
#spe_IF_out <- here("processed-data", "spot_deconvo", "05-shared_utilities", "IF", "spe.h5ad")
spe_HE_out <- here("processed-data", "spot_deconvo", "shared_utilities", "HE", "spe.h5ad")

#sample_IF_out <- here("processed-data", "spot_deconvo", "05-shared_utilities", "IF", "sample_ids.txt")
sample_HE_out <- here("processed-data", "spot_deconvo", "shared_utilities", "HE", "sample_ids.txt")

marker_object_out <- here("processed-data", "spot_deconvo", "shared_utilities", paste0("marker_stats_", cell_group, ".rds"))

if (cell_group == "broad") {
  cell_type_var <- "broad.type"
} else {
  cell_type_var <- "cell.type"
}

#  Make sure output directories exist
# dir.create(dirname(spe_IF_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(spe_HE_out), recursive = TRUE, showWarnings = FALSE)

###############################################################################
#  Functions
###############################################################################

write_anndata <- function(sce, out_path) {
  invisible(
    basiliskRun(
      fun = function(sce, filename) {
        library("zellkonverter")
        library("reticulate")
        
        # Convert SCE to AnnData:
        adata <- SCE2AnnData(sce)
        
        #  Write AnnData object to disk
        adata$write(filename = filename)
        
        return()
      },
      env = zellkonverterAnnDataEnv(),
      sce = sce,
      filename = out_path
    )
  )
}

###############################################################################
#   Main
###############################################################################

#   Load objects
load(sce_in, verbose = TRUE)

# spe_IF <- readRDS(spe_IF_in)

load(spe_HE_in, verbose = TRUE)
spe_HE <- spe
rm(spe)
gc()

print(paste0("Running script at ", cell_group, "-resolution."))

#-------------------------------------------------------------------------------
#   Drop appropriate cells
#-------------------------------------------------------------------------------

# if (cell_group == "layer") {
#   #   Drop poor/suspicious clusters and subclusters
#   keep <- !is.na(sce$layer_level) & (sce$cellType_hc != "drop")
# } else {
#   #   Drop suspicious subclusters
#   keep <- sce$cellType_hc != "drop"
# }
# 
# print("Distribution of cells to drop (FALSE) vs. keep (TRUE):")
# table(keep)
# sce <- sce[, keep]

#-------------------------------------------------------------------------------
#   Add cell counts to HE spatial object
#-------------------------------------------------------------------------------
# HE_id_path <- here("processed-data", "spot_deconvo", "HE_ID_table.csv")
# HE_counts_path <- here("processed-data", "rerun_spaceranger", "{sample_id}", "outs", "spatial", "tissue_spot_counts.csv")

# id_table <- read.csv(HE_id_path)
# 
# spe_HE$count <- NA
# for (sample_id in unique(spe_HE$sample_id)) {
#   #   Correctly determine the path for the cell counts for this sample, then
#   #   read in
#   long_id <- id_table[match(sample_id, id_table$short_id), "long_id"]
#   this_path <- sub("{sample_id}", long_id, HE_counts_path, fixed = TRUE)
#   cell_counts <- read.csv(this_path)
#   
#   #   All spots in the object should have counts
#   stopifnot(
#     all(
#       colnames(spe_nonIF[, spe_nonIF$sample_id == sample_id]) %in%
#         cell_counts$barcode
#     )
#   )
#   
#   #   Line up the rows of 'cell_counts' with the sample-subsetted SPE object
#   cell_counts <- cell_counts[
#     match(
#       colnames(spe_nonIF[, spe_nonIF$sample_id == sample_id]),
#       cell_counts$barcode
#     ),
#   ]
#   
#   #   Add this sample's counts to the SPE object
#   spe_nonIF$count[spe_nonIF$sample_id == sample_id] <- cell_counts$Nmask_dark_blue
# }
# 
# #   Ensure counts were read in for all spots in the object
# if (any(is.na(spe_nonIF$count))) {
#   stop("Did not find cell counts for all non-IF spots.")
# }

#-------------------------------------------------------------------------------
#   Convert snRNA-seq and spatial R objects to AnnData python objects
#-------------------------------------------------------------------------------

#   zellkonverter doesn't know how to convert the 'spatialCoords' slot. We'd
#   ultimately like the spatialCoords in the .obsm['spatial'] slot of the
#   resulting AnnDatas, which corresponds to reducedDims(spe)$spatial in R
#reducedDims(spe_IF)$spatial <- spatialCoords(spe_IF)
reducedDims(spe_HE)$spatial <- spatialCoords(spe_HE)

#   Use Ensembl gene IDs for rownames (not gene symbol)
rownames(sce) <- rowData(sce)$gene_id

#   Save a copy of the filtered + slightly modified sce as an R object, and
#   convert all objects to Anndatas
saveRDS(sce, sce_r_out)

print("Converting objects to AnnDatas...")
write_anndata(sce, sce_out)

#   Spatial objects are the same between broad and layer-level resolutions, and
#   need only be saved once
if (cell_group == "broad") {
  # write_anndata(spe_IF, spe_IF_out)
  write_anndata(spe_HE, spe_HE_out)
  
  #   Write sample names to text files
  # writeLines(unique(spe_IF$sample_id), con = sample_IF_out)
  # writeLines(unique(spe_nonIF$sample_id), con = sample_nonIF_out)
}
gc()

#-------------------------------------------------------------------------------
#   Rank marker genes
#-------------------------------------------------------------------------------

#   We won't consider genes that aren't in both the spatial objects
# keep <- (rowData(sce)$gene_id %in% rowData(spe_IF)$gene_id) &
#   (rowData(sce)$gene_id %in% rowData(spe_HE)$gene_id)
keep = (rowData(sce)$gene_id %in% rowData(spe_HE)$gene_id)
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

saveRDS(marker_stats, marker_object_out)

session_info()
