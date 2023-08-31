setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(spatialLIBD)
    library(here)
    library(edgeR)
    library(scuttle)
    library(scater)
    library(scran)
    library(dplyr)
    library(PCAtools)
    library(sessioninfo)
    library(gridExtra)
    library(ggforce)
    library(ggspavis)
})
csv<-read.csv('hpc_allSamples_graphst_k16.csv')
# Remove last dash and the number after it using gsub
csv$barcode <- gsub("-[0-9]+$", "", csv$spot_id)
csv$sample_id<-factor(csv$captureArea)
levels(csv$sample_id)[33:36]<-gsub(levels(csv$sample_id)[33:36],
                                   pattern='V12F14-051',
                                   replacement='Br2720')
csv$key<-paste0(csv$barcode,'_',as.character(csv$sample_id))

#now rearrange csv so it matches the colnames of spe
csv<-csv[match(colnames(spe),csv$key),]
#and add clusters to spe
colnames(csv)
colData(spe)$graphst<-csv$cluster_lamb_0_1_and_1
spe$graphst<-factor(spe$graphst)

##broad groupings
spe_pseudo$broad.class<-factor(ifelse(spe_pseudo$graphst %in% c(2,3,4,7,9,11,14,15),'Neuron',
                        ifelse(spe_pseudo$graphst %in% c(6,8),'Neuropil/WM',
                        ifelse(spe_pseudo$graphst %in% c(5,12,13),'WM',
                        ifelse(spe_pseudo$graphst %in% c(16),'artifacts?',
                        ifelse(spe_pseudo$graphst %in% c(5,12,13),'WM',
                        'Vascular/CSF'))))))
