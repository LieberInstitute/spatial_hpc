setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(spatialLIBD)
    library(here)
    library(scater)
    library(scran)
    library(nnSVG)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(sessioninfo)
})

load(file = here::here("processed-data", "06_Clustering", "spe_modify.Rdata"))
spe <- speB
# Load BayesSpace clusters onto spe object
spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data", "06_Clustering", "BayesSpace"),
    prefix = ""
)

# Create vector of samples for nnSVG on whole tissue
brains <- as.character(unique(spe$brnum))
samples <- unique(spe$sample_id)
# Run nnSVG once per sample whole tissue and store lists of top SVGs
res_list <- as.list(rep(NA, length(brains)))

names(res_list) <- brains

for (s in seq_along(brains)) {

    # select sample_id
    # ix <- spe$brnum == brains[s]
    ix <- spe$sample_id == samples[s]
    spe_sub <- spe[, ix]

    # run nnSVG filtering for mitochondrial gene and low-expressed genes
    spe_sub <- filter_genes(spe_sub)

    # re-calculate logcounts after filtering
    spe_sub <- logNormCounts(spe_sub)

    # run nnSVG
    set.seed(12345)
    spe_sub <- nnSVG(spe_sub, n_threads = 10)

    # store whole tissue results
    res_list[[s]] <- rowData(spe_sub)

    save(res_list[[s]], file = here::here("processed-data","07_Feature_selection", "nnSVG", paste(spe$sample_id, ".RData")))
}

# save whole tissue nnSVG results
save(res_list, file = here::here("processed-data","07_Feature_selection", "nnSVG", "nnSVG.RData"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
