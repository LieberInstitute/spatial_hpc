library(RcppML)
library(pheatmap)
library(SingleCellExperiment)
load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))
###Run NMF
x<-RcppML::nmf(assay(sce,'logcounts'),
    k=100,
    tol = 1e-06,
    maxit = 1000,
    verbose = T,
    L1 = 0.1,
    seed = 1135,
    mask_zeros = FALSE,
    diag = TRUE,
    nonneg = TRUE
)

##add pats to sce
loads<-t(x@h)
colData(sce)<-cbind(colData(sce),loads)

save(sce,file=here::here('snRNAseq_hpc','processed-data','sce','sce_nmf_final.rda'))
save(x,file=here::here('snRNAseq_hpc','processed-data','nmf','nmf_final.rda'))
