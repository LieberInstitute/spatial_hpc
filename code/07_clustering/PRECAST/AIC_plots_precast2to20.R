
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

##write loop to get all the ICs
kvals=2:20
###load final PRECAST clusters
aicList<-list()
for(k in kvals){
print('loading object ',k)
load(file = here("processed-data", "06_clustering", "PRECAST",
                 paste0("nnSVG_2000_",k,"_HE.Rdata")))
PRECASTObj <- selectModel(PRECASTObj,criteria='AIC')
aicList[k]<-PRECASTObj@resList$icMat[2]
print(aicList[k])
rm(PRECASTObj)
gc()
}

