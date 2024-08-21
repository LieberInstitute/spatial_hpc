library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("here")
library("batchelor")
library(dplyr)

load("snRNAseq_hpc/python_analysis/processed-data/sce_postqc-loose.Rdata")
#rownames(sce)<-uniquifyFeatureNames(rowData(sce)$gene_id,rowData(sce)$gene_name)

##make sure these are factors
sce$brnum<-factor(sce$brnum)
sce$round<-factor(sce$round)
sce$sorted<-factor(sce$sorted)
# feature selection
set.seed(800)
sce <- devianceFeatureSelection(sce,
                                assay = "counts", fam = "poisson",
                                sorted = T,batch=sce$brnum)
hdg<-rownames(counts(sce))[1:2000]
res<-sce[rownames(counts(sce)) %in% hdg,]

#null residuals
set.seed(913)
res <- nullResiduals(res, fam = "poisson", type = "pearson")

#PCA
set.seed(915)
res <- scater::runPCA(res,ntop=2000,
                      exprs_values='poisson_pearson_residuals',
                      ncomponents=100,
                      BSPARAM=BiocSingular::IrlbaParam())

reducedDim(sce,'PCA')<-reducedDim(res,'PCA')
save(sce, file="snRNAseq_hpc/python_analysis/processed-data/sce_postqc-loose.Rdata")
