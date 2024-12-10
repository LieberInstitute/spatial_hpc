setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc")
# check if Eriks code for nmf is reproducible (identical results)
# source code: https://github.com/LieberInstitute/spatial_hpc/blob/main/snRNAseq_hpc/code/nmf/nmf.R
suppressPackageStartupMessages({
library(RcppML)
library(pheatmap)
library(SingleCellExperiment)
library(Matrix)
})

load(file="snRNAseq_hpc/processed-data/sce/sce_final.rda")
###Run NMF
newx<-RcppML::nmf(assay(sce,'logcounts'),
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
saveRDS(newx, "snRNAseq_hpc/processed-data/revision/nmf_redo.rda")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
