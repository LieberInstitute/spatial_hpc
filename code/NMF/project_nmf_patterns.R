library(CoGAPS)
library(RcppML)
library(SpatialExperiment)
load(file=here::here(NMF_PATH_HERE))
load(file=here::here(SPE_PATH_HERE)

proj<-projectR(
    data=data,
    loadings=x@w,
    full = FALSE,
    family = "gaussianff",
    bootstrapPval = FALSE,
    bootIter = 1000
)

set.seed(1029)
i<-intersect(rownames(spe),rownames(x@w))
loadings<-x@w
loadings<-loadings[rownames(loadings) %in% i,]
spe2<-spe[rownames(spe) %in% i,]
loadings<-loadings[match(rownames(spe2),rownames(loadings)),]

proj<-project(loadings,logcounts(spe2),L2=0.00005)
