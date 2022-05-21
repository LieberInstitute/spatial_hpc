library("spatialLIBD")
library("here")

## Load the spe object

load(file=here::here("processed-data", "02_build_spe", "spe_basic.Rdata", verbose = TRUE))
#load("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/02_build_spe/spe_basic.Rdata")
lobstr::obj_size(spe) / 1024^3
# 4.856404 B
imgData(spe) <- imgData(spe)[!imgData(spe)$image_id %in% c("hires", "detected", "aligned"),]
lobstr::obj_size(spe) / 1024^3
# 2.132348 B
save(spe, file = here::here("processed-data","02_build_spe","spe_small.Rdata"))
