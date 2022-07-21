setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(SpatialExperiment)
    library(spatialLIBD)
    library(BayesSpace)
    library(tidySingleCellExperiment)
})

load(file = here::here("processed-data", "05_Batch_correction", "spe_harmony.Rdata"))
spe <- cluster_import(
    spe,
    cluster_dir = here::here("processed-data", "06_Clustering", "BayesSpace"),
    prefix = ""
)

table(spe$BayesSpace_harmony_k11_nrep10000, spe$brnum)

#     Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325 Br8492 Br8667
# 1    1683   2064   1743   1094   1636   3034   2499   1127   2005
# 2     781   1957   1270    724    750   1166   2431    335   1244
# 3      45   1179     61     65   1560    245    879   1401   3202
# 4    4601   1741   5379   1582   1682   3071   2959    109   1006
# 5    1693   2238   1673   1269   2754   3602   2536   1206   2870
# 6    1277   2388   1078    888   1279   2471   1663    688   1316
# 7    1499   2604    806   1412    726    945   1863    841   2820
# 8    1983   3473   1137   1124   1530   1920   3492    703   1409
# 9      68    817     66    226   1304    448    647   1203    986
# 10    282    590    479    150    739    549    716    568    468
# 11    175    389    311     85    271   1189    730    191    507

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
save(spe, file = here::here("processed-data", "06_Clustering", "spe_modify.Rdata"))
