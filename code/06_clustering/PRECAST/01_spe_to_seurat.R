setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages({
library("here")
library("SpatialExperiment")
library("Seurat")
library("SeuratData")
})

load(file = here::here("processed-data", "05_preprocess_batchCorrection", "spe_norm.rda"))
source(file = here::here("code", "06_clustering", "BayesSpace", "preprocess_harmony", "offset_check.R"))
x<-spe
#x = offset_check(spe)
#
### check
# df <- cbind.data.frame(colData(x), spatialCoords(x))
# ggplot(df, aes(x = row, y = col, color = sample_id)) +
# geom_point(size = 1) +
# coord_fixed() +
# guides(color = guide_legend(override.aes = list(size = 3))) +
# theme_bw()
#  
#x$spot_id = x$key
#colData(x)$row<-colData(x)$array_row
#colData(x)$col<-colData(x)$array_col

sue <- CreateSeuratObject(
      counts=as.matrix(counts(x)),
      meta.data=data.frame(colData(x)),
      project="HPC")

seuList = list()

samples = unique(x$sample_id)
table(x$brnum, x$sample_id)

for (i in seq_along(samples)){
seuList[[i]] = subset(x=sue, subset = sample_id == samples[i])
}

save(seuList, file = here("processed-data", "06_clustering", "PRECAST", "seuList_counts.Rdata"))

sue <- CreateSeuratObject(
      counts=as.matrix(logcounts(x)),
      meta.data=data.frame(colData(x)),
      project="HPC")

seuList = list()

samples = unique(x$sample_id)
table(x$brnum, x$sample_id)

for (i in seq_along(samples)){
seuList[[i]] = subset(x=sue, subset = sample_id == samples[i])
}

save(seuList, file = here("processed-data", "06_clustering", "PRECAST", "seuList_logcounts.Rdata"))


source(file = here::here("code", "06_clustering", "BayesSpace", "preprocess_harmony", "offset_check.R"))
x<-spe
x = offset_check(spe)

## check
 df <- cbind.data.frame(colData(x), spatialCoords(x))
 ggplot(df, aes(x = row, y = col, color = sample_id)) +
 geom_point(size = 1) +
 coord_fixed() +
 guides(color = guide_legend(override.aes = list(size = 3))) +
 theme_bw()
  
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

save(seuList, file = here("processed-data", "06_clustering", "PRECAST", "seuList_counts_brnum.Rdata"))
