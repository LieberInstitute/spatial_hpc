
setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("lobstr"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("ggspavis"))

load(file = here::here("processed-data", "02_build_spe", "spe_raw.Rdata"))

dim(spe_raw)
# [1]  36601 159744

sampleID <- unique(spe_raw$sample_id)
# [1] "V10B01-086_C1" "V10B01-086_D1" "V11U08-081_C1" "V11U08-081_D1"
# [5] "V11L05-333_A1" "V11L05-333_B1" "V11L05-333_C1" "V11L05-333_D1"
# [9] "V10B01-085_A1" "V10B01-085_B1" "V10B01-085_C1" "V10B01-085_D1"
# [13] "V10B01-086_A1" "V10B01-086_B1" "V11L05-335_A1" "V11L05-335_B1"
# [17] "V11L05-335_C1" "V11U08-084_A1" "V11U08-084_B1" "V11U08-084_C1"
# [21] "V11U08-084_D1" "V11A20-297_A1" "V11A20-297_B1" "V11A20-297_C1"
# [25] "V11A20-297_D1" "V11L05-335_D1" "V11U08-081_A1" "V11U08-081_B1"
# [29] "V11L05-336_A1" "V11L05-336_B1" "V11L05-336_C1" "V11L05-336_D1"

# Rotations by sample:
angle_list <- c(0, 180, 0, 0,
                0, 0, 0, 0,
                90, 180, 180, 180, 
                0, 0, 0, 0, 
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0)
source(file = here::here("code", "02_build_spe", "transform_spe.R"))

for (i in seq_along(angle_list)) {
    id <- sampleID[i]
    x <- trans_geom(spe_raw, sample_id = id, degrees = angle_list[i])

    if (i == 1) {
        spe <- x
    } else {
        spe <- cbind(spe, x)
    }

    spe <- rotateImg(spe, sample_id = id, image_id = TRUE, angle_list[i])
   }


pdf(file = here::here("plots", "02_build_spe", "ReferenceMapping1.pdf"), h = 10, w = 20)
for (i in seq_along(sampleID)) {
    id <- sampleID[i]
    x_raw = spe_raw[, (colData(spe_raw)$in_tissue & colData(spe_raw)$sample_id == id)]
    x_trans = spe[, (colData(spe)$in_tissue & colData(spe)$sample_id == id)]
    p1 <- plotVisium(x_raw, spots = TRUE, y_reverse = TRUE)
    p2 <- plotVisium(x_trans, spots = TRUE, y_reverse = TRUE)
    grid.arrange(p1, p2, ncol = 2)
    
    df_raw <- cbind.data.frame(colData(x_raw), spatialCoords(x_raw))
    df_trans <- cbind.data.frame(colData(x_trans), spatialCoords(x_trans))
    p3 = ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + 
         geom_point(size = 1) + coord_fixed() + theme_bw()
    p4 = ggplot(df_trans, aes(x = array_row, y = array_col, color = sample_id)) + 
        geom_point(size = 1) + coord_fixed() + theme_bw() 
    
    grid.arrange(p3, p4, ncol = 2)
}
dev.off()

table(spe$sample_id,spe$position)

#                      BL   BR    N   TL   TR
# V10B01-085_A1 4992    0    0    0    0    0
# V10B01-085_B1 4992    0    0    0    0    0
# V10B01-085_C1 4992    0    0    0    0    0
# V10B01-085_D1 4992    0    0    0    0    0
# V10B01-086_A1    0    0    0    0 4992    0
# V10B01-086_B1    0 4992    0    0    0    0
# V10B01-086_C1    0    0    0    0    0 4992
# V10B01-086_D1    0    0 4992    0    0    0
# V11A20-297_A1    0 4992    0    0    0    0
# V11A20-297_B1    0    0 4992    0    0    0
# V11A20-297_C1    0    0    0    0 4992    0
# V11A20-297_D1    0    0    0    0    0 4992
# V11L05-333_A1    0    0    0    0 4992    0
# V11L05-333_B1    0    0    0    0    0 4992
# V11L05-333_C1    0 4992    0    0    0    0
# V11L05-333_D1    0    0 4992    0    0    0
# V11L05-335_A1    0    0 4992    0    0    0
# V11L05-335_B1    0    0    0    0    0 4992
# V11L05-335_C1    0    0    0    0 4992    0
# V11L05-335_D1    0    0    0 4992    0    0
# V11L05-336_A1    0    0 4992    0    0    0
# V11L05-336_B1    0 4992    0    0    0    0
# V11L05-336_C1    0    0    0    0    0 4992
# V11L05-336_D1    0    0    0    0 4992    0
# V11U08-081_A1    0    0    0    0 4992    0
# V11U08-081_B1    0    0    0    0    0 4992
# V11U08-081_C1    0    0    0    0 4992    0
# V11U08-081_D1    0    0    0    0    0 4992
# V11U08-084_A1    0    0    0    0 4992    0
# V11U08-084_B1    0    0    0    0    0 4992
# V11U08-084_C1    0 4992    0    0    0    0
# V11U08-084_D1    0    0 4992    0    0    0

unique(spe$sample_id[which(spe$position == "N")])
spe$position[which(spe$sample_id == "V11L05-335_D1")] <- "BR"

unique(spe$sample_id[which(spe$position == "")])
spe$position[which(spe$sample_id == "V10B01-085_A1")] <- "TR"
spe$position[which(spe$sample_id == "V10B01-085_B1")] <- "TL"
spe$position[which(spe$sample_id == "V10B01-085_C1")] <- "BR"
spe$position[which(spe$sample_id == "V10B01-085_D1")] <- "BL"
spe$position <- factor(spe$position, levels = c("TL", "TR", "BR", "BL"))

