
setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("cowplot"))

###############################################################################
#  Visually check quality of markers
###############################################################################
source(here("code","spot_deconvo","shared_utilities","plottingfunctions.R"))
source(here("code", "spot_deconvo", "shared_utilities", "shared_function.R"))

Dr <- here("processed-data","spot_deconvo","shared_utilities")

cell_group = "broad"
cell_type_var = "broad.class"
name = "_class"

cell_group = "layer" 
cell_type_var = "cell.class"
name = "_celltype_class1_noHATAGABAAmy"

n_markers_per_type <- 50
cell_type_nrow <- 4

plot_dir <- here("plots", "spot_deconvo", "shared_utilities", cell_group)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

#   Load objects
#sce = readRDS(here(Dr,"sce.rds"))
#sce = readRDS(here(Dr,"sce_class.rds"))
sce = readRDS(here(Dr,"sce_class_noHATAGABA.rds"))
spe = readRDS(here(Dr,"spe.rds"))
#spg = readRDS(here(Dr,"spg.rds"))
marker_stats = readRDS(here(Dr,paste0("marker_stats_",cell_group,name,".rds")))
cell_types = unique(marker_stats$cellType.target)

#   Visually show how markers look for each cell type
plot_list <- lapply(
  cell_types,
  function(ct) {
    genes <- marker_stats |>
      filter(
        rank_ratio <= n_markers_per_type,
        cellType.target == ct,
        ratio > 1
      ) |>
      pull(gene)
    my_plotExpression(
      sce, genes,
      ct = cell_type_var,
      #fill_colors = metadata(sce)[[colors_col]],
      title = paste("Top", length(genes), "for", ct),
      marker_stats = marker_stats |>
        filter(
          rank_ratio <= n_markers_per_type,
          cellType.target == ct,
          ratio > 1
        )
    )
  }
)

#   Write a multi-page PDF with violin plots for each cell group and all markers
pdf(file.path(plot_dir, paste0("marker_gene_violin",cell_group,name,".pdf")),width = 35, height = 35)
print(plot_list)
dev.off()

#   Plot mean ratio against log fold-change for all genes, split by target cell
#   type and colored by whether each gene will be used as a marker
p <- marker_stats |>
  mutate(
    Marker = case_when(
      rank_ratio <= n_markers_per_type ~ paste0(
        "Top-", n_markers_per_type, " Marker"
      ),
      TRUE ~ "Not Marker"
    )
  ) |>
  ggplot(aes(ratio, std.logFC, color = Marker)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  facet_wrap(~cellType.target, scales = "free_x", nrow = cell_type_nrow) +
  labs(x = "Mean Ratio") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  guides(col = guide_legend(override.aes = list(size = 2)))

pdf(file.path(plot_dir, paste0("mean_ratio_vs_1vall",cell_group,name,".pdf")), width = 10)
print(p)
dev.off()

#   Plot mean-ratio distibution by group (cell type or layer label)
boxplot_mean_ratio(n_markers_per_type, "mean_ratio_boxplot")

#   Get Ensembl ID for classical markers
# spT <- spg
classical_markers <- c("PPFIA2", "AMPH", "FNDC1", "GFRA1", "KRT17", "C5orf63", "GAD2", "MIF", "FABP7", "MAN1A2", "SFRP2", "MOBP", "MAG", "MTURN", "PHLDB1", "ACTA2", "TTR")
clusters <-c("GCL", "CA2-4", "CA1", "SUB", "SUB.RHP", "RHP", "GABA", "SL_SR", "ML", "SR_SLM", "SLM_WM", "WM", "WM1", "WM2", "WM3", "Vascular", "Choroid")

rownames(sce) <- rowData(sce)$gene_name
spT <- spe[, which(spe$brnum == "Br3942")] #use this temporarily until we get SPG data
stopifnot(all(classical_markers %in% rowData(sce)$gene_name))
classical_markers_ens <- rownames(sce)[match(classical_markers, rowData(sce)$gene_name)]
stopifnot(all(classical_markers_ens %in% rowData(spT)$gene_name))
rownames(spT) <- rowData(spT)$gene_name
#-------------------------------------------------------------------------------
#   First, plot the spatial expression of classical markers as reference
#-------------------------------------------------------------------------------

plot_list <- list()
i <- 1

for (j in 1:length(classical_markers)) {
  for (sample_id in unique(spT$sample_id)) {
    #   Determine the title for this subplot
      title <- paste0(classical_markers[j], ": marker for cluster ", clusters[j], "\n(", sample_id, ")")
    #   Produce the ggplot object (grid version)
    plot_list[[i]] <- spot_plot(
      spT,
      sample_id = sample_id, var_name = classical_markers_ens[j],
      include_legend = TRUE, is_discrete = FALSE, title = title,
      assayname = "counts", minCount = 0
    )
    
    i <- i + 1
  }
}

write_spot_plots(
  plot_list = plot_list, n_col = length(unique(spT$sample_id)),
  plot_dir = plot_dir, file_prefix = "marker_spatial_sparsity_reference",
  include_individual = FALSE
)

#-------------------------------------------------------------------------------
#   For IF, show a grid of plots summarizing how sparsely marker genes
#   for each cell type are expressed spatially. Repeat these plots for different
#   numbers of markers per cell type: 15, 25, 50
#-------------------------------------------------------------------------------
#rownames(spT)=rowData(spT)$gene_name
rownames(spT)=rowData(spT)$gene_id # for deconvobuddies

for (n_markers in c(15, 25)) {
  plot_list <- list()
  i <- 1
  
  #   Plot proportion of markers having nonzero expression for each cell type
  for (ct in cell_types) {
    #   Get markers for this cell type
    markers <- marker_stats |>
      filter(cellType.target == ct,
        rank_ratio <= n_markers,
        ratio > 1) |>
      pull(gene)
    # markers = rownames(marksList[[ct]])
    # keep = (markers %in% rowData(spT)$gene_name)
    # markers = (markers[keep])[1:n_markers]
    for (sample_id in unique(spT$sample_id)) {
        spT_small <- spT[markers, spT$sample_id == sample_id]
      
      #   For each spot, compute proportion of marker genes with nonzero expression
      spT_small$prop_nonzero_marker <- colMeans(assays(spT_small)$counts > 0)
      
      plot_list[[i]] <- spot_plot(
        spT_small,
        sample_id = sample_id,
        var_name = "prop_nonzero_marker", include_legend = TRUE,
        is_discrete = FALSE, minCount = 0,
        title = paste0("Prop. markers w/ nonzero exp:\n", ct, " (", sample_id, ")")
      )
      
      i <- i + 1
    }
  }
  n_sample <- length(unique(spT$sample_id))
  #n_rows <- length(unique(marker_stats$cellType.target))
  
  write_spot_plots(
    plot_list = plot_list, n_col = n_sample, plot_dir = plot_dir,
    file_prefix = paste0("marker_spatial_sparsity_n", n_markers),
    include_individual = TRUE
  )
}

#-------------------------------------------------------------------------------
#   For the manuscript, make a 4x4 plot for IF layer-level: 4 columns (samples)
#   by 4 rows: expression of PCP4 then that of 'Excit_L5' for different numbers
#   of markers: 15,25,50
#-------------------------------------------------------------------------------
rownames(spT)=rowData(spT)$gene_name
if (cell_group == "layer") {
  plot_list <- list()
  max_list <- list()
  i <- 1
  
  #   Plot expression of PCP4 for every sample
  for (sample_id in unique(spT$sample_id)) {
    #   Produce the ggplot object
    plot_list[[i]] <- spot_plot(
      spT,
      sample_id = sample_id,
      var_name = classical_markers_ens[classical_markers == "PPFIA2"],
      include_legend = TRUE, is_discrete = FALSE, title = "PPFIA2 Counts",
      assayname = "counts", minCount = 0
    )
    
    max_list[[i]] <- 0
    
    i <- i + 1
  }
  
  for (n_markers in c(15, 25, 50)) {
    #   Get markers for this cell type
    markers <- marker_stats |>
      filter(
        cellType.target == "GC",
        rank_ratio <= n_markers,
        ratio > 1
      ) |>
      pull(gene)
    
    for (sample_id in unique(spT$sample_id)) {
      spe_small <- spT[markers, spT$sample_id == sample_id]
      
      #   For each spot, compute proportion of marker genes with nonzero
      #   expression
      spe_small$prop_nonzero_marker <- colMeans(
        assays(spe_small)$counts > 0
      )
      
      max_list[[i]] <- max(spe_small$prop_nonzero_marker)
      
      plot_list[[i]] <- spot_plot(
        spe_small,
        sample_id = sample_id,
        var_name = "prop_nonzero_marker", include_legend = TRUE,
        is_discrete = FALSE, minCount = 0,
        title = paste0(
          "Prop. GC markers w/ >0 counts (n = ", n_markers, ")"
        )
      )
      
      i <- i + 1
    }
  }
  
  max_mat <- matrix(
    unlist(max_list),
    ncol = length(unique(spe$sample_id)), byrow = TRUE
  )
  
  #   Now loop back through the plot list (which will be displayed in 2D)
  #   and overwrite the scale to go as high as the largest value in the
  #   column. This allows for easy comparison between number of markers for
  #   the same samples
  for (i_col in 1:length(unique(spe$sample_id))) {
    for (i_row in 2:4) {
      index <- (i_row - 1) * length(unique(spe$sample_id)) + i_col
      upper_limit <- max(max_mat[, i_col])
      
      plot_list[[index]] <- overwrite_scale(
        plot_list[[index]],
        upper_limit = upper_limit, min_count = 0
      )
    }
  }
  
  write_spot_plots(
    plot_list = plot_list, n_col = n_sample, plot_dir = plot_dir,
    file_prefix = "sparsity_figure", include_individual = FALSE
  )
}

session_info()