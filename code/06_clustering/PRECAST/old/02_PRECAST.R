setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages({
library(dplyr)
library(purrr)
library(Seurat)
library(SpatialExperiment)
library(PRECAST)
library(tictoc)
library("here")
})

load(file = here("processed-data", "06_clustering", "PRECAST", "PRECASTObj_17.Rdata"))

resList <- PRECASTObj@resList
PRECASTObj <- selectModel(PRECASTObj)
# What is ARI statistics?
# ari_precast <- mclust::adjustedRandIndex(PRECASTObj@resList$cluster[[1]], PRECASTObj@seulist[[1]]$layer_guess_reordered)

seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")
seuInt

# Edit by ADR, changed color palette from cols_cluster with Polychome (Maddy
# already started this change)

# cols_cluster <- c("#FD7446", "#709AE1", "#31A354", "#9EDAE5", "#DE9ED6", "#BCBD22", "#CE6DBD", "#DADAEB",
#                   "yellow", "#FF9896", "#91D1C2", "#C7E9C0", "#6B6ECF", "#7B4173", "red")
cols <- Polychrome::palette36.colors(16)
names(cols) <- sort(unique(seuInt$cluster))

p12 <- SpaPlot(seuInt, batch = NULL, item = "cluster", point_size = 2, cols = cols)

pdf(file = here::here("plots", "06_clustering", "PRECAST", "K16(1).pdf"), width = 21, height = 20)
for (i in c(1:9)){
  print(p12[[i]])
}
dev.off()

p12

seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")
tmp_colData <- colData(x) |> as.data.frame() |>
  tibble::rownames_to_column(var = "spot_id") |>
  left_join(seuInt@meta.data |> select(batch, cluster) |>
              tibble::rownames_to_column(var = "spot_id"))

fnl_spe_sub <- spe_sub
tmp_colData <- colData(fnl_spe_sub) |> as.data.frame() |>
  left_join(seuInt@meta.data |> select(batch, cluster)|> tibble::rownames_to_column(var = "spot_id"))

colData(fnl_spe_sub)$batch <- tmp_colData$batch
colData(fnl_spe_sub)$cluster <- tmp_colData$cluster

library(ggspavis)
# library(RColorBrewer)
nb.cols <- 18
mycolors <- colorRampPalette(unname(palette.colors(36, "Polychrome 36")))(nb.cols)

plotVisium(fnl_spe_sub, fill = "cluster", palette = mycolors, alpha = 1)
ggsave("hbc_plots.png")
Footer
