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

k_vals<-15:25

#preList<-list()
#for(k in k_vals){
#i=k-14
#print(paste0('loading PRECASTObj: k=',k))
#load(file = here("processed-data", "06_clustering", "PRECAST",
#                 paste0("PRECASTObj_nnSVG_2000_",k,"brnum.Rdata")))
#preList[[i]]<-PRECASTObj@resList
#rm(PRECASTObj)
#gc()
#                 }
#


aicList<-list()
print('starting loop!')
for(k in k_vals){
i=k-14
print(paste0('loading PRECASTObj: k=',k))
load(file = here("processed-data", "06_clustering", "PRECAST",
                 paste0("PRECASTObj_nnSVG_2000_",k,"_counts.Rdata")))

# Leaving this bit from Maddy's script, not sure I need it.
resList <- PRECASTObj@resList
PRECASTObj <- SelectModel(PRECASTObj,criteria='AIC')
aicList[i]<-PRECASTObj@resList$icMat[2]
print(aicList[i])
print('Seurat conversion')
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")
seuInt

dim(spe)

dim(seuInt@meta.data)

# Find number of barcodes that are different from spe

length(setdiff(colnames(spe), rownames(seuInt@meta.data)))

# Create temp df for the clusters from PRECAST k_15

k_tab <- data.frame(
    cluster = seuInt@meta.data$cluster,
    barcodes = rownames(seuInt@meta.data)
)

rownames(k_tab) <- rownames(seuInt@meta.data)

#indices <- which(rownames(k_tab) %in% setdiff(rownames(k_tab), colnames(spe)))
#k_tab<-k_tab[-indices,]

# Create df of NAs to fill missing barcodes & in an reorder barcodes

#emptyNaDF <- data.frame(matrix(NA,nrow = length(setdiff(colnames(spe), rownames(seuInt@meta.data))), ncol = 2))

#rownames(emptyNaDF) <- setdiff(colnames(spe), rownames(seuInt@meta.data))

#colnames(emptyNaDF) <- c("cluster", "barcodes")

#k_tab <- rbind(k_tab, emptyNaDF)

print('matching')
k_tab <- k_tab[order(match(rownames(k_tab),colnames(spe))), ]

stopifnot(rownames(k_tab) == colnames(spe))

k_label<-paste0("PRECAST_sample_k",k)

print('adding to SPE')
spe[[k_label]] <- factor(k_tab$cluster)
print(paste0('finished loop iteration',k))
rm(list = setdiff(ls(), c("spe", "k_vals","aicList",'hk_genes')))

}

print('loop over, saving SPE')
save(spe,file=here::here('processed-data','02_build_spe','spe_precast.rda'))
