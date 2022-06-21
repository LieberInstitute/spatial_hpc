# cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/

library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("jaffelab")
library("tidyverse")
library("here")
library("sessioninfo")
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("ggspavis"))
suppressPackageStartupMessages(library("gridExtra"))

#### Compute QC metrics ####
load(here("processed-data", "04_QC", "spe_QC.Rdata"), verbose = TRUE)

spebr6522 = spe[,which(spe$brnum == "Br6522")]

spebr6522 <- scuttle::addPerCellQC(
  spebr6522,
  subsets = list(Mito = which(seqnames(spebr6522) == "chrM")),
  BPPARAM = BiocParallel::MulticoreParam(4)
)

#### Check for low quality spots ####

## High mito
spebr6522$high_mito_id <- isOutlier(spebr6522$subsets_Mito_percent, nmads = 1, type = "higher", batch = spebr6522$sample_id)
spebr6522$high_mito_br <- isOutlier(spebr6522$subsets_Mito_percent, nmads = 1, type = "higher", batch = spebr6522$brnum)

pdf(here("plots","04_QC", "QC_Br6522.pdf"), width = 21, height = 10)
plotColData(spebr6522, x = "brnum", y = "subsets_Mito_percent", colour_by = "high_mito_br") +
  ggtitle("Mito Precent Br6522, nmad=5") + 
  facet_wrap(~ spebr6522$sample_id, scales = "free_x", nrow = 1)

samples = unique(spebr6522$sample_id)
# QC plot of tissue spots discarded
for (i in 1:4){
  p = vis_clus(
    spe = spebr6522,
    sampleid = samples[i],
    clustervar = "high_mito_br",
    colors = c("FALSE" = "grey90", "TRUE" = "red"),
    point_size = 2,
    ... = paste0("_Br6522")
  )
  
  p1 = plotVisium(spebr6522[,which(spebr6522$sample_id == samples[i])], spots = FALSE)
  
  grid.arrange(p,p1,nrow = 1)
}

for (i in 1:4){
  p = vis_clus(
    spe = spebr6522,
    sampleid = samples[i],
    clustervar = "high_mito_id",
    colors = c("FALSE" = "grey90", "TRUE" = "red"),
    point_size = 2,
    ... = paste0("_Br6522_by capture area")
  )
  
  p1 = plotVisium(spebr6522[,which(spebr6522$sample_id == samples[i])], spots = FALSE)
  
  grid.arrange(p,p1,nrow = 1)
}

dev.off()


vis_grid_clus(
  spe = spe,
  clustervar = "high_mito_br",
  pdf = here::here("plots","04_QC","QC_discard_mito_by_brain.pdf"),
  spatial = FALSE,
  sort_clust = FALSE,
  colors = c("FALSE" = "grey90", "TRUE" = "red"),
  point_size = 2
)

vis_grid_clus(
  spe = spe,
  clustervar = "high_mito_id",
  pdf = here::here("plots","04_QC","QC_discard_mito_by_capture_area.pdf"),
  spatial = FALSE,
  sort_clust = FALSE,
  colors = c("FALSE" = "grey90", "TRUE" = "red"),
  point_size = 2
)
