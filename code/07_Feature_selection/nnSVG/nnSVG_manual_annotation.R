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

# Create vector of samples for nnSVG on whole tissue
brains <- as.character(unique(spe$brnum))
samples <- unique(spe$sample_id)
# Run nnSVG once per sample whole tissue and store lists of top SVGs
res_list <- as.list(rep(NA, length(samples)))

names(res_list) <- samples

for (s in seq_along(samples)) {
  
  # select sample_id
  # ix <- spe$brnum == brains[s]
  ix <- spe$sample_id == samples[s]
  spe_sub <- spe[, ix]
  spe_sub <- spe_sub[, which(spe_sub$ManualAnnotation == "GCL")]
  
  # run nnSVG filtering for mitochondrial gene and low-expressed genes
  spe_sub <- filter_genes(spe_sub)
  
  # re-calculate logcounts after filtering
  spe_sub <- logNormCounts(spe_sub)
  
  # run nnSVG
  set.seed(12345)
  message('running nnSVG')
  spe_sub <- nnSVG(spe_sub, n_threads = 10)
  
  # store whole tissue results
  message('saving data')
  res_list[[s]] <- rowData(spe_sub)
  temp = rowData(spe_sub)
  save(temp, file = here::here("processed-data","07_Feature_selection", "nnSVG_manual_annotation", paste0(samples[s], ".Rdata")))
}

# save whole tissue nnSVG results
save(res_list, file = here::here("processed-data","07_Feature_selection", "nnSVG_manual_annotation", "nnSVG.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
