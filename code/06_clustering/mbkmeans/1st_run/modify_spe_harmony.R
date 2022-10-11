setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(SpatialExperiment)
    library(spatialLIBD)
    library(BayesSpace)
    library(tidySingleCellExperiment)
})

load(file = here::here("processed-data", "05_preprocess_batchCorrection", "spe_harmony.Rdata"))
spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data", "06_Clustering", "BayesSpace", "1st_run"),
    prefix = ""
)


## change colnames(spe) and rownames(spatialCoords(spe)) to spe$key
colnames(spe) <- spe$key
rownames(spatialCoords(spe)) <- spe$key
## now arrange() to get the correct order
spe <- arrange(spe, slide, array)

samples <- unique(spe$sample_id)
angle_list <- c(90, 180, 180, 180, 0, 0, 270, 0, 0, 0, 0, 0, 0, 0, 0, 180, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
source(file = here::here("code", "pilot_data_checks", "transform_spe.R"))

for (i in seq_along(angle_list)) {
    id <- samples[i]
    x <- trans_geom(spe, sample_id = id, degrees = angle_list[i])
    if (i == 1) {
        speB <- x
    } else {
        speB <- cbind(speB, x)
    }
    speB <- rotateImg(speB, sample_id = id, image_id = TRUE, angle_list[i])
}

speB$position <- factor(speB$position, levels = c("TL", "TR", "BL", "BR"))
speB <- arrange(speB, brnum, position)
spe = speB
save(spe, file = here::here("processed-data", "06_Clustering", "BayesSpace", "1st_run", "spe_modify.Rdata"))
