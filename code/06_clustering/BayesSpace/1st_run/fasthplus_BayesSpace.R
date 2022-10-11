setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library(fasthplus)
library(SpatialExperiment)
library(here)
library(sessioninfo)
library(spatialLIBD)

k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# load spe object
load(file = here::here("processed-data", "06_Clustering", "BayesSpace", "1st_run", "spe_modify.Rdata"), verbose = TRUE)

dim(spe)
# [1] 30359 135640

# hpb estimate. t = pre-bootstrap sample size, D = reduced dimensions matrix, L = cluster labels, r = number of bootstrap iterations
dim(reducedDims(spe)$HARMONY)
# [1] 135640     50

find_t <- function(L, proportion = 0.05) {
  initial_t <- floor(length(L) * proportion)
  smallest_cluster_size <- min(table(L))
  n_labels <- length(unique(L))
  ifelse(smallest_cluster_size > (initial_t / n_labels), initial_t, smallest_cluster_size * n_labels)
}

message(paste0("Find fasthplus for clusters ",k, " - ") , Sys.time())
initial_t <- find_t(L = colData(spe)[[paste0("BayesSpace_harmony_k", k, "_nrep10000")]], proportion = 0.01)

set.seed(20220216)
fasthplus <- hpb(D = reducedDims(spe)$HARMONY, L = colData(spe)[[paste0("BayesSpace_harmony_k", k, "_nrep10000")]], t = initial_t, r = 30)
results <- data.frame(k = k, fasthplus = fasthplus)
write.table(results, file = here::here("processed-data", "06_Clustering", "fasthplus_results1.csv"), append = TRUE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
Footer
