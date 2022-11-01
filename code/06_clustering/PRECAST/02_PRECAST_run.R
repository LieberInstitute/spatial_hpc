setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(Seurat)
  library(SpatialExperiment)
  library(PRECAST)
  library(tictoc)
})

load(file = here::here("processed-data", "06_clustering", "PRECAST", "seuList.Rdata"))

preobj = CreatePRECASTObject(seuList = seuList)
preobj@seulist

PRECASTObj <- AddAdjList(preobj, platform = "Visium")
## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the
## information in the algorithm.
PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, coreNum = 1, maxIter = 30, verbose = TRUE)

tic()
PRECASTObj <- PRECAST(PRECASTObj, K = 15)
toc()

save(PRECASTObj, file = here("processed-data", "06_clustering", "PRECAST", "PRECASTObj.Rdata"))
