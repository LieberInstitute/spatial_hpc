setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages({
  library("dplyr")
  library("purrr")
  library("Seurat")
  library("here")
  library("sessioninfo")
  library("SpatialExperiment")
  library("PRECAST")
  library("tictoc")
})

load(file = here::here("processed-data", "06_clustering", "PRECAST", "seuList_allSamples.Rdata"))

preobj = CreatePRECASTObject(seuList = seuList)
preobj@seulist

PRECASTObj <- AddAdjList(preobj, platform = "Visium")
## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the
## information in the algorithm.
PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, coreNum = 8, maxIter = 30, verbose = TRUE)

K <- as.numeric(Sys.getenv("SGE_TASK_ID"))

tic()
PRECASTObj <- PRECAST(PRECASTObj, K = K)
toc()

save(PRECASTObj, file = here("processed-data", "06_clustering", "PRECAST", paste0("allSamples_PRECASTObj_",K,".Rdata")))
