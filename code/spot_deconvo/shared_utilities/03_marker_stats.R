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
