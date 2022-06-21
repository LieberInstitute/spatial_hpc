# cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("ggspavis"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("here"))

load(here("processed-data", "04_QC", "spe_QC.Rdata"), verbose = TRUE)

samples <- unique(colData(spe)[,c("sample_id","brnum")])
rownames(samples) <- NULL
pdf(here("plots","04_QC", "QC_discard_brain.pdf"), width = 21, height = 10)

# QC plot of tissue spots discarded
for (i in 1:32){
p = vis_clus(
  spe = spe,
  sampleid = samples$sample_id[i],
  clustervar = "discard_auto_br",
  colors = c("FALSE" = "grey90", "TRUE" = "red"),
  point_size = 2,
  ... = paste0("_",samples$brnum[i])
)

p1 = plotVisium(spe[,which(spe$sample_id == samples$sample_id[i])], spots = FALSE)

grid.arrange(p,p1,nrow = 1)
}

dev.off()

pdf(here("plots","04_QC", "QC_discard_capture_area.pdf"), width = 21, height = 10)
# QC plot of tissue spots discarded

for (i in 1:32){
  p = vis_clus(
    spe = spe,
    sampleid = samples$sample_id[i],
    clustervar = "discard_auto_id",
    colors = c("FALSE" = "grey90", "TRUE" = "red"),
    point_size = 2,
    ... = paste0("_",samples$brnum[i])
  )
  
  p1 = plotVisium(spe[,which(spe$sample_id == samples$sample_id[i])], spots = FALSE)
  
  grid.arrange(p,p1,nrow = 1)
}

dev.off()

vis_grid_clus(
  spe = spe,
  clustervar = "discard_auto_id",
  pdf = here::here("plots","04_QC","QC_discard_all_capture_area.pdf"),
  spatial = FALSE,
  sort_clust = FALSE,
  colors = c("FALSE" = "grey90", "TRUE" = "red"),
  point_size = 2
)

vis_grid_clus(
  spe = spe,
  clustervar = "discard_auto_br",
  pdf = here::here("plots","04_QC","QC_discard_all_brain.pdf"),
  spatial = FALSE,
  sort_clust = FALSE,
  colors = c("FALSE" = "grey90", "TRUE" = "red"),
  point_size = 2
)
