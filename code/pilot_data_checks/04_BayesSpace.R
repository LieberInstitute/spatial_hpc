
setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
    library("here")
    library("spatialLIBD")
    library("SpatialExperiment")
    library("ggplot2")
    library("patchwork")
    library("scater")
    library("harmony")
    library("BayesSpace")
    library("scran")
})


load(file = here::here("processed-data", "pilot_data_checks", "spe_harmony.Rdata"))
# OFFSET
# The spatial locations for each sample need to be offset so that spots of different samples are not neighbors.
# summary(spatialData(speH)$array_row)
# summary(spatialData(speH)$array_col)
auto_offset_row <- as.numeric(factor(unique(speH$sample_id))) * 100
names(auto_offset_row) <- unique(speH$sample_id)
speH$row <- spatialData(speH)$array_row + auto_offset_row[speH$sample_id]
speH$col <- spatialData(speH)$array_col
# summary(colData(speH)$row)
# summary(colData(speH)$col)

# pdf(file=here::here("plots", "pilot_data_checks", "hpc_BayesSpace_OffsetCheck.pdf"))
# clusterPlot(speH, "subject", color = NA) + #make sure no overlap between samples
#   labs(fill = "Subject", title = "Offset check")
# dev.off()
#
# Sys.time()
# speB = spatialCluster(speH, use.dimred = "HARMONY", q = 7, nrep = 10000) #use HARMONY
# Sys.time()
#
#
# pdf(file=here::here("plots", "pilot_data_checks", "hpc_BayesSpace_clusterPlot_10k.pdf"))
# clusterPlot(speB, color = NA) + #plot clusters
#   labs(title = "BayesSpace joint clustering")
# dev.off()
#
# save(speB, file=here::here("processed-data","pilot_data_checks", "spe_bayesSpace_10k.Rdata"))

Sys.time()
speB <- spatialCluster(speH, use.dimred = "HARMONY", q = 7, nrep = 50000) # use HARMONY
Sys.time()

pdf(file = here::here("plots", "pilot_data_checks", "hpc_BayesSpace_clusterPlot_50k.pdf"))
clusterPlot(speB, color = NA) + # plot clusters
    labs(title = "BayesSpace joint clustering")
dev.off()

save(speB, file = here::here("processed-data", "pilot_data_checks", "spe_bayesSpace_50k.Rdata"))

load(file = here::here("processed-data", "pilot_data_checks", "spe_bayesSpace_50k.Rdata"))

p_list <- vis_grid_clus(
    speB[, speB$subject %in% c("Br6423-O")],
    "spatial.cluster",
    spatial = FALSE,
    return_plots = TRUE,
    sort_clust = FALSE,
    point_size = 2,
    "Br6423-O"
)

pdf(file = here::here("plots", "pilot_data_checks", "hpc_BayesSpace_50k_visgrid_Br6423-O.pdf"), h = 8, w = 12)
cowplot::plot_grid(plotlist = p_list, ncol = 2)
dev.off()

p_list1 <- vis_grid_clus(
    speB[, speB$subject %in% c("Br2743-Y")],
    "spatial.cluster",
    spatial = FALSE,
    return_plots = TRUE,
    sort_clust = FALSE,
    point_size = 2,
    "Br2743-Y"
)

p_list2 <- vis_grid_clus(
    speB[, speB$subject %in% c("Br6432-R")],
    "spatial.cluster",
    spatial = FALSE,
    return_plots = TRUE,
    sort_clust = FALSE,
    point_size = 2,
    "Br6432-R"
)

pdf(file = here::here("plots", "pilot_data_checks", "hpc_BayesSpace_50k_visgrid.pdf"), h = 8, w = 12)
cowplot::plot_grid(plotlist = c(p_list1, p_list2), ncol = 2)
dev.off()

## vis_grid for cluster.init
p_list <- vis_grid_clus(
    speB[, speB$subject %in% c("Br6423-O")],
    "cluster.init",
    spatial = FALSE,
    return_plots = TRUE,
    sort_clust = FALSE,
    point_size = 2,
    "Br6423-O"
)

pdf(file = here::here("plots", "pilot_data_checks", "hpc_BayesSpace_50k_init_Br6423-O.pdf"), h = 8, w = 12)
cowplot::plot_grid(plotlist = p_list, ncol = 2)
dev.off()

p_list1 <- vis_grid_clus(
    speB[, speB$subject %in% c("Br2743-Y")],
    "cluster.init",
    spatial = FALSE,
    return_plots = TRUE,
    sort_clust = FALSE,
    point_size = 2,
    "Br2743-Y"
)

p_list2 <- vis_grid_clus(
    speB[, speB$subject %in% c("Br6432-R")],
    "cluster.init",
    spatial = FALSE,
    return_plots = TRUE,
    sort_clust = FALSE,
    point_size = 2,
    "Br6432-R"
)

pdf(file = here::here("plots", "pilot_data_checks", "hpc_BayesSpace_50k_init.pdf"), h = 8, w = 12)
cowplot::plot_grid(plotlist = c(p_list1, p_list2), ncol = 2)
dev.off()
