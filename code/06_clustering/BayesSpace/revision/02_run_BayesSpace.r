setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(BayesSpace)
})
# vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/BayesSpace/inst/doc/BayesSpace.html
### needs log transformed data
### needs PCs
set.seed(123)

load("processed-data/06_clustering/BayesSpace/revision/spe_bayes.Rdata")
spe$row <- spe$array_row
spe$col <- spe$array_col

cat("\n",format(Sys.time()),"spatialPreprocess\n")
spe <- spatialPreprocess(spe, platform="Visium", skip.PCA=TRUE, log.normalize=FALSE)

cat(format(Sys.time()),"spatialCluster\n")
spe <- spatialCluster(spe, q=18, use.dimred="MNN", d=20, #start with using 20 MNN dims
               init.method="kmeans", #initial clustering method options = mclust or kmeans
               model="t", #Using t-distributed errors in likelihood/bayesian models is a different way to obtain robust methods, as the t-distribution has heavier tails than the normal.
               nrep=10000 #default nrep is 50k, vignette recommends min 10k so will start there
               )
table(colData(spe)[,c("cluster","spatial.cluster")])

cat("\n\n",format(Sys.time()),"save results\n")
save(spe, file="processed-data/06_clustering/BayesSpace/revision/spe_bayes_k18-kmeans-10k.Rdata")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
