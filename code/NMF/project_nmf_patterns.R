####snRNA-seq NMF pattern projection to Visium
library(RcppML)
library(SpatialExperiment)
library(here)

##load nmf patterns
load(file=here::here('snRNAseq_hpc','processed-data','NMF','nmf_final.rda'))

##load spe
load(file=here::here('processed-data','06_clustering',
                     'PRECAST','spe_precast_HE_domain.rda'))

##projection
set.seed(1029)
i<-intersect(rownames(spe),rownames(x@w))
loadings<-x@w
loadings<-loadings[rownames(loadings) %in% i,]
spe2<-spe[rownames(spe) %in% i,]
loadings<-loadings[match(rownames(spe2),rownames(loadings)),]

proj<-project(loadings,logcounts(spe2),L2=0.00005)

###rescale
proj<-apply(proj,1,function(x){x/sum(x)})

###cbind coldata and projections
colData(spe)<-cbind(colData(spe),proj)

###save data
save(spe,file=here::here('processed-data','NMF',
                     'spe_nmf_final.rda'))
