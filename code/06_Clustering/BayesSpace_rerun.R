###############################
# Set up SGE array job to run k=5 to k = 20
# Found in BayesSpaces.sh shell script line -t 2-15
setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
  library(here)
  library(sessioninfo)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(BayesSpace)
  library(ggplot2)
  library(Polychrome)
})

# Load SPE
load(file = here::here("processed-data", "06_Clustering", "spe_modify.Rdata"))

# Choose k
k <- 17

load(file = here("processed-data", "06_Clustering", "mbkmeans.Rdata"))
spe$kmeans <- km_res[[13]]$Clusters

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

# The spatial locations for each sample need to be offset so that spots of different samples are not neighbors.
# summary(colData(spe)$array_row)
# summary(colData(spe)$array_col)
auto_offset_row <- as.numeric(factor(unique(spe$sample_id))) * 100
names(auto_offset_row) <- unique(spe$sample_id)
colData(spe)$row <- colData(spe)$array_row + auto_offset_row[spe$sample_id]
colData(spe)$col <- colData(spe)$array_col

# Run BayesSpace
message("Running spatialCluster()")
Sys.time()
set.seed(12345)
spe <- spatialCluster(spe, use.dimred = "HARMONY", q = k, platform = "Visium", init.method = "kmeans", save.chain = TRUE, nrep = 10000)
Sys.time()

save(spe, file = here("processed-data", "06_Clustering", "BayesSpace_rerun_k17_mbkmeans.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
