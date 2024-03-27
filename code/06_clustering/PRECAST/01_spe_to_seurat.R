setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages({
library("here")
library("SpatialExperiment")
library("Seurat")
library("SeuratData")
})

load(file = here::here("processed-data", "05_preprocess_batchCorrection", "spe_norm_final.rda"))
spe<-spe[,spe$brnum %in% levels(spe$brnum)[1:10]]
spe$brnum<-droplevels(spe$brnum)

source(file = here::here("code", "06_clustering", "PRECAST", "offset_check1.R"))
x = offset_check(spe)

### check
pdf(file=here::here('plots','06_clustering','PRECAST','offset_check_HE.pdf'))
 df <- cbind.data.frame(colData(x), spatialCoords(x))
 ggplot(df, aes(x = row, y = col, color = sample_id)) +
 geom_point(size = 1) +
 coord_fixed() +
 guides(color = guide_legend(override.aes = list(size = 3))) +
 theme_bw()
dev.off()
 
x$spot_id = x$key
 
sue <- CreateSeuratObject(
  counts=as.matrix(counts(x)),
  meta.data=data.frame(colData(x)),
  project="HPC")
 
seuList = list()
 
brains = unique(x$brnum)
table(x$brnum, x$sample_id)
 
for (i in seq_along(brains)){
  seuList[[i]] = subset(x=sue, subset = brnum == brains[i])
 }
 
save(seuList, file = here("processed-data", "06_clustering", "PRECAST", "seuList_allSamples.Rdata"))