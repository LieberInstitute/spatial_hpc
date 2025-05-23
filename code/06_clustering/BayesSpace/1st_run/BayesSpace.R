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

# Create directory for BayesSpace plots
dir_plots <- here::here("plots", "06_Clustering", "BayesSpace", "1st_run")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
load(file = here::here("processed-data", "05_preprocess_batchCorrection", "spe_harmony.Rdata"))

# Choose k
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

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
spe <- spatialCluster(spe, use.dimred = "HARMONY", q = k, platform = "Visium", save.chain = TRUE, nrep = 10000)
Sys.time()

nrep <- 10000
spe$BayesSpace_temp <- spe$spatial.cluster
BayesSpace_name <- paste0("BayesSpace_harmony_k", k, "_nrep", nrep)
colnames(colData(spe))[ncol(colData(spe))] <- BayesSpace_name

cluster_export(
    spe,
    BayesSpace_name,
    cluster_dir = here::here("processed-data", "06_Clustering", "BayesSpace", "1st_run")
)

## Visualize BayesSpace results
cols <- Polychrome::palette36.colors(k)
names(cols) <- sort(unique(spe$spatial.cluster))

vis_grid_clus(
    spe = spe,
    clustervar = paste0("BayesSpace_harmony_k", k, "_nrep", nrep),
    pdf_file = here("plots", "06_Clustering", "BayesSpace", "1st_run", paste0("vis_grid_clus_BayesSpace_k", k, "_nrep", nrep, ".pdf")),
    sort_clust = FALSE,
    colors = cols,
    spatial = FALSE,
    point_size = 1,
    sample_order = unique(spe$sample_id)
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
