setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("lobstr"))
suppressPackageStartupMessages(library("sessioninfo"))

load(file = here::here("processed-data", "02_build_spe", "spe_transform.Rdata"))

dim(spe)
# [1]  36601 219648
lobstr::obj_size(spe)
# 8.04 GB

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)
length(no_expr)
# [1] 5118
length(no_expr) / nrow(spe) * 100
# [1] 13.98322
spe <- spe[-no_expr, ]
dim(spe)
# [1]  31483 219648

## Now drop the spots outside the tissue
spe <- spe[, colData(spe)$in_tissue]
dim(spe)
# [1]  31483 191140

## Remove spots without counts
if (any(colSums(counts(spe)) == 0)) {
  message("removing spots without counts for spe")
  spe <- spe[, -which(colSums(counts(spe)) == 0)]
  dim(spe)
}

# [1]  31483 191136


save(spe, file = here::here("processed-data", "02_build_spe", "spe_basic.Rdata"))

# session_info()
