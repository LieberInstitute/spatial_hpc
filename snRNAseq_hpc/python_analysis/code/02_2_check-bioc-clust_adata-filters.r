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

# remaining steps performed in Rstudio
##MNN batch correction
#set.seed(1788)
#mnn<-batchelor::reducedMNN(reducedDim(sce,'PCA'),batch=sce$brnum,k=50)
#reducedDim(sce,'MNN')<-mnn$corrected
##set.seed(1788)
#mnn<-batchelor::reducedMNN(reducedDim(sce,'MNN'),batch=sce$round,k=50)
#reducedDim(sce,'MNN')<-mnn$corrected

##umap
#sce <- runUMAP(sce, dimred = "MNN")

#sce = addPerCellQC(sce)
#sce$lg10.sum = log10(sce$sum)

# comparison to erik's final clusters also performed in rstudio
#sce_keep <- sce
#load(here("snRNAseq_hpc",'processed-data','sce','sce_final.rda'))
#cmp.df = left_join(as.data.frame(colData(sce_keep)), as.data.frame(colData(sce)[,1:31]),
#                   by=c("Sample","Barcode","brnum","round","sorted"="sort"))
#filt.df = filter(cmp.df, !is.na(cell.type))
#nrow(filt.df)

#sce_keep$erik_final = FALSE
#colData(sce_keep)[sce_keep$key %in% filt.df$key.x, "erik_final"] = TRUE
#sce_keep$cell.type="none"
#colData(sce_keep)[sce_keep$erik_final==T,"cell.type"] = as.character(filt.df$cell.type)

#sce <- sce_keep
#save(sce, file="snRNAseq_hpc/python_analysis/processed-data/sce_postqc-strict.Rdata")
