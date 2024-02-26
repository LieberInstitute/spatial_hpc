setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(Seurat)
    library(SpatialExperiment)
    library(PRECAST)
    library(spatialLIBD)
    library(ggplot2)
    library(gridExtra)
    library("here")
})
print('Loading spe!')
load(file="/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/05_preprocess_batchCorrection/spe_norm.rda")
spe<-spe[,spe$brnum %in% levels(spe$brnum)[1:10]]


print(paste0('loading PRECASTObj'))
load(file = here("processed-data", "06_clustering", "PRECAST",
                 "PRECASTObj_nnSVG_18_final.Rdata"))

#Format for addition to spe
resList <- PRECASTObj@resList
PRECASTObj <- selectModel(PRECASTObj)

print('Seurat conversion')
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")
seuInt

dim(spe)

dim(seuInt@meta.data)

# Find number of barcodes that are different from spe

length(setdiff(colnames(spe), rownames(seuInt@meta.data)))

# Create temp df for the clusters 

k_tab <- data.frame(
    cluster = seuInt@meta.data$cluster,
    barcodes = rownames(seuInt@meta.data)
)

rownames(k_tab) <- rownames(seuInt@meta.data)

print('matching')
k_tab <- k_tab[order(match(rownames(k_tab),colnames(spe))), ]

stopifnot(rownames(k_tab) == colnames(spe))

k_label<-paste0("PRECAST_k18")

print('adding to SPE')
spe[[k_label]] <- factor(k_tab$cluster)

save(spe,file=here::here('processed-data','06_clustering','PRECAST','spe_precast_HE.rda'))
