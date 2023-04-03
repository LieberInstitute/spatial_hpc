# set up arrayjob to run k=11 to k = 20
# don't use spatial preprocess. in order to do this you have to reset metadata
# increase nrep for spatialCluster??
setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library("BayesSpace")
library("RColorBrewer")
library("ggplot2")
library("gridExtra")
})

temp = "OSCApreprocess_harmony_captureArea" #"spatialPreprocess_harmony", "OSCApreprocess_harmony_captureArea", "OSCApreprocess_harmony_brain"
## load data
load(file = here::here("processed-data", "05_preprocess_batchCorrection", paste0(temp ,"_spe_allSamples.Rdata")), verbose = TRUE)
dim(spe)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))


source(file = here::here("code", "06_clustering", "BayesSpace", "preprocess_harmony", "offset_check.R"))
spe = offset_check(spe)

## check
 df <- cbind.data.frame(colData(x), spatialCoords(x))
 ggplot(df, aes(x = row, y = col, color = sample_id)) +
 geom_point(size = 1) +
 coord_fixed() +
 guides(color = guide_legend(override.aes = list(size = 3))) +
 theme_bw()

### BayesSpace on Batch Corrected
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

message("Running spatialCluster()")
Sys.time()
set.seed(20220201)
spe <- spatialCluster(spe, use.dimred = "HARMONY", q = k,nrep=10000)
Sys.time()

#spe <- spatialCluster(spe, use.dimred = "HARMONY", q = k, nrep = 10000)

spe$bayesSpace_temp <- spe$spatial.cluster
bayesSpace_name <- paste0("bayesSpace_captureArea_", k)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
  spe,
  bayesSpace_name,
  cluster_dir = here::here("processed-data", "06_clustering", "BayesSpace", "preprocess_harmony", bayesSpace_name)
)

sample_ids <- unique(colData(spe)$sample_id)
cols <- Polychrome::palette36.colors(k)
names(cols) <- sort(unique(spe$spatial.cluster))
brains = unique(spe$brnum)
clustV = bayesSpace_name
  
pdf(file = here::here("plots", "06_clustering", "BayesSpace", "preprocess_harmony", paste0(bayesSpace_name, ".pdf")), width = 21, height = 20)

for (i in seq_along(brains)){
  speb <- spe[, (colData(spe_raw)$brnum == brains[i])]
  samples <- unique(speb$sample_id)
  samples
  
  if (length(samples) == 1){
    p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = clustV, colors = cols, point_size = 2, ... = paste0("_", brains[i]))
  } else if (length(samples) == 2){
    p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = clustV, colors = cols, point_size = 2, ... = paste0("_", brains[i]))
    p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = clustV, colors = cols, point_size = 2, ... = paste0("_", brains[i]))
    grid.arrange(p1, p2, nrow = 2)
  } else if (length(samples) == 3){
    p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = clustV, colors = cols, point_size = 2, ... = paste0("_", brains[i]))
    p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = clustV, colors = cols, point_size = 2, ... = paste0("_", brains[i]))
    p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = clustV, colors = cols, point_size = 2, ... = paste0("_", brains[i]))
    grid.arrange(p1, p2, p3, nrow = 2)
  } else if (length(samples) == 4){
    p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = clustV, colors = cols, point_size = 2, ... = paste0("_", brains[i]))
    p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = clustV, colors = cols, point_size = 2, ... = paste0("_", brains[i]))
    p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = clustV, colors = cols, point_size = 2, ... = paste0("_", brains[i]))
    p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = clustV, colors = cols, point_size = 2, ... = paste0("_", brains[i]))
    grid.arrange(p1, p2, p3, p4, nrow = 2)
  } else if (length(samples) == 5){
    p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = clustV, colors = cols, point_size = 2, ... = paste0("_", brains[i]))
    p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = clustV, colors = cols, point_size = 2, ... = paste0("_", brains[i]))
    p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = clustV, colors = cols, point_size = 2, ... = paste0("_", brains[i]))
    p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = clustV, colors = cols, point_size = 2, ... = paste0("_", brains[i]))
    p5 <- vis_clus(spe = speb, sampleid = samples[5], clustervar = clustV, colors = cols, point_size = 2, ... = paste0("_", brains[i]))
    grid.arrange(p1, p2, p3, p4, p5, nrow = 2)}
}
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()