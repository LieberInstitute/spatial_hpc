library(RcppML)
library(singlet)
library(here)
library(sessioninfo)
setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')

load(here::here('snRNAseq_hpc','processed-data','NMF','logcounts.rda'))
cvnmf<-cross_validate_nmf(
    mtx,
    ranks=c(5,10,50,100,125,150,200),
    n_replicates = 3,
    tol = 1e-03,
    maxit = 100,
    verbose = 3,
    L1 = 0.1,
    L2 = 0,
    threads = 0,
    test_density = 0.2
)

save(cvnmf,file=here::here('snRNAseq_hpc','processed-data','NMF','cvnmf_3reps.rda'))
