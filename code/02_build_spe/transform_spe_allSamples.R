#cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/

suppressPackageStartupMessages({
library("here")
library("SpatialExperiment")
library("spatialLIBD")
library("rtracklayer")
library("lobstr")
library("sessioninfo")
library("ggplot2")
library("gridExtra")
library("ggspavis")
library("tidySingleCellExperiment")
library("dplyr")
library("tidyverse")
})


load(file = here::here("processed-data", "02_build_spe", "spe_raw_allSamples.Rdata"))
spe_raw = spe
dim(spe_raw)
# [1]  36601 179712

spe_raw$sample_id <- paste0(spe_raw$slide,"_",spe_raw$array)
sampleID <- unique(spe_raw$sample_id)
sampleID
#  [1] "V10B01-086_C1" "V10B01-086_D1" "V11U08-081_C1" "V11U08-081_D1"
#  [5] "V11L05-333_A1" "V11L05-333_B1" "V11L05-333_C1" "V11L05-333_D1"
#  [9] "V10B01-085_A1" "V10B01-085_B1" "V10B01-085_C1" "V10B01-085_D1"
#  [13] "V10B01-086_A1" "V10B01-086_B1" "V11L05-335_A1" "V11L05-335_B1"
#  [17] "V11L05-335_C1" "V11U08-084_A1" "V11U08-084_B1" "V11U08-084_C1"
#  [21] "V11U08-084_D1" "V11A20-297_A1" "V11A20-297_B1" "V11A20-297_C1"
#  [25] "V11A20-297_D1" "V11L05-335_D1" "V11U08-081_A1" "V11U08-081_B1"
#  [29] "V11L05-336_A1" "V11L05-336_B1" "V11L05-336_C1" "V11L05-336_D1"
#  [33] "V12F14-051_A1" "V12F14-051_B1" "V12F14-051_C1" "V12F14-051_D1"

colnames(spe_raw) <- spe_raw$key
rownames(spatialCoords(spe_raw)) <- spe_raw$key

## arrange by brain and array (A1, B1, C1, D1)
spe_raw <- arrange(spe_raw, brnum, array)

samples <- unique(spe_raw$sample_id)
samples
# [1] "V10B01-086_C1" "V11U08-081_C1" "V10B01-086_D1" "V11U08-081_D1"
# [5] "V11L05-333_A1" "V11L05-333_B1" "V11L05-333_C1" "V11L05-333_D1"
# [9] "V10B01-085_A1" "V10B01-085_B1" "V10B01-085_C1" "V10B01-085_D1"
# [13] "V10B01-086_A1" "V10B01-086_B1" "V11L05-335_A1" "V11L05-335_B1"
# [17] "V11L05-335_C1" "V11U08-084_A1" "V11U08-084_B1" "V11U08-084_C1"
# [21] "V11U08-084_D1" "V11A20-297_A1" "V11A20-297_B1" "V11A20-297_C1"
# [25] "V11A20-297_D1" "V11L05-335_D1" "V11U08-081_A1" "V11U08-081_B1"
# [29] "V11L05-336_A1" "V11L05-336_B1" "V11L05-336_C1" "V11L05-336_D1"

brains = unique(spe_raw$brnum)
pdf(here("plots", "02_build_spe", "prearrangedAllSamples.pdf"), width = 8, height = 10)

for (i in seq_along(brains)){
    speb <- spe_raw[, which(spe_raw$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    samples
    
    if (length(samples) == 1){
        p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = FALSE)
    } else if (length(samples) == 2){
        p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = FALSE)
        p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = FALSE)
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = FALSE)
        p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = FALSE)
        p3 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = FALSE)
        grid.arrange(p1, p2, p3, nrow = 2)
    } else if (length(samples) == 4){
        p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = FALSE)
        p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = FALSE)
        p3 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = FALSE)
        p4 <- plotVisium(speb[, which(speb$sample_id == samples[4])], spots = FALSE)
        grid.arrange(p1, p2, p3, p4, nrow = 2)
    } else if (length(samples) == 5){
        p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = FALSE)
        p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = FALSE)
        p3 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = FALSE)
        p4 <- plotVisium(speb[, which(speb$sample_id == samples[4])], spots = FALSE)
        p5 <- plotVisium(speb[, which(speb$sample_id == samples[5])], spots = FALSE)
        grid.arrange(p1, p2, p3, p4, p5, nrow = 2)}
}

dev.off()

table(spe_raw$sample_id,spe_raw$brnum)

# Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325 Br8492 Br8667
# V10B01-085_A1      0      0   4992      0      0      0      0      0      0
# V10B01-085_B1      0      0   4992      0      0      0      0      0      0
# V10B01-085_C1      0      0   4992      0      0      0      0      0      0
# V10B01-085_D1      0      0   4992      0      0      0      0      0      0
# V10B01-086_A1      0      0      0   4992      0      0      0      0      0
# V10B01-086_B1      0      0      0   4992      0      0      0      0      0
# V10B01-086_C1   4992      0      0      0      0      0      0      0      0
# V10B01-086_D1   4992      0      0      0      0      0      0      0      0
# V11A20-297_A1      0      0      0      0      0      0   4992      0      0
# V11A20-297_B1      0      0      0      0      0      0   4992      0      0
# V11A20-297_C1      0      0      0      0      0      0   4992      0      0
# V11A20-297_D1      0      0      0      0      0      0   4992      0      0
# V11L05-333_A1      0   4992      0      0      0      0      0      0      0
# V11L05-333_B1      0   4992      0      0      0      0      0      0      0
# V11L05-333_C1      0   4992      0      0      0      0      0      0      0
# V11L05-333_D1      0   4992      0      0      0      0      0      0      0
# V11L05-335_A1      0      0      0      0   4992      0      0      0      0
# V11L05-335_B1      0      0      0      0   4992      0      0      0      0
# V11L05-335_C1      0      0      0      0   4992      0      0      0      0
# V11L05-335_D1      0      0      0      0      0      0   4992      0      0
# V11L05-336_A1      0      0      0      0      0      0      0      0   4992
# V11L05-336_B1      0      0      0      0      0      0      0      0   4992
# V11L05-336_C1      0      0      0      0      0      0      0      0   4992
# V11L05-336_D1      0      0      0      0      0      0      0      0   4992
# V11U08-081_A1      0      0      0      0      0      0      0   4992      0
# V11U08-081_B1      0      0      0      0      0      0      0   4992      0
# V11U08-081_C1   4992      0      0      0      0      0      0      0      0
# V11U08-081_D1   4992      0      0      0      0      0      0      0      0
# V11U08-084_A1      0      0      0      0      0   4992      0      0      0
# V11U08-084_B1      0      0      0      0      0   4992      0      0      0
# V11U08-084_C1      0      0      0      0      0   4992      0      0      0
# V11U08-084_D1      0      0      0      0      0   4992      0      0      0
# V12F14-051_A1      0      0      0      0      0      0      0      0      0
# V12F14-051_B1      0      0      0      0      0      0      0      0      0
# V12F14-051_C1      0      0      0      0      0      0      0      0      0
# V12F14-051_D1      0      0      0      0      0      0      0      0      0
# 
# Br2720
# V10B01-085_A1      0
# V10B01-085_B1      0
# V10B01-085_C1      0
# V10B01-085_D1      0
# V10B01-086_A1      0
# V10B01-086_B1      0
# V10B01-086_C1      0
# V10B01-086_D1      0
# V11A20-297_A1      0
# V11A20-297_B1      0
# V11A20-297_C1      0
# V11A20-297_D1      0
# V11L05-333_A1      0
# V11L05-333_B1      0
# V11L05-333_C1      0
# V11L05-333_D1      0
# V11L05-335_A1      0
# V11L05-335_B1      0
# V11L05-335_C1      0
# V11L05-335_D1      0
# V11L05-336_A1      0
# V11L05-336_B1      0
# V11L05-336_C1      0
# V11L05-336_D1      0
# V11U08-081_A1      0
# V11U08-081_B1      0
# V11U08-081_C1      0
# V11U08-081_D1      0
# V11U08-084_A1      0
# V11U08-084_B1      0
# V11U08-084_C1      0
# V11U08-084_D1      0
# V12F14-051_A1   4992
# V12F14-051_B1   4992
# V12F14-051_C1   4992
# V12F14-051_D1   4992

samples = unique(spe_raw$sample_id)
samples
# Br2743 [1] "V10B01-086_C1" "V11U08-081_C1" "V10B01-086_D1" "V11U08-081_D1"
# Br3942 [5] "V11L05-333_A1" "V11L05-333_B1" "V11L05-333_C1" "V11L05-333_D1"
# Br6423 [9] "V10B01-085_A1" "V10B01-085_B1" "V10B01-085_C1" "V10B01-085_D1"
# Br6432 [13] "V10B01-086_A1" "V10B01-086_B1" 
# Br6471 [15] "V11L05-335_A1" "V11L05-335_B1" "V11L05-335_C1" 
# Br6522 [18] "V11U08-084_A1" "V11U08-084_B1" "V11U08-084_C1" "V11U08-084_D1" 
# Br8325 [22] "V11A20-297_A1" "V11A20-297_B1" "V11A20-297_C1" "V11A20-297_D1"  "V11L05-335_D1" 
# Br8492 [27] "V11U08-081_A1" "V11U08-081_B1"
# Br8667 [29] "V11L05-336_A1" "V11L05-336_B1" "V11L05-336_C1" "V11L05-336_D1"
# Br2720 [33] "V12F14-051_A1" "V12F14-051_B1" "V12F14-051_C1" "V12F14-051_D1"

position_list <- c("TR", "BL", "TL", "BR",
                   "TL", "TR", "BL", "BR",
                   "TR", "TL", "BR", "BL",
                   "TL", "BL",
                   "BL", "TR", "TL",
                   "TL", "TR", "BL", "BR",
                   "BL", "BR", "TL", "TR", "TR",
                   "TL", "BL",
                   "TL", "TR", "BL", "BR",
                   "BL","BR","TL","TR")

for (i in seq_along(samples)){
    spe_raw$position[colData(spe_raw)$sample_id == samples[i]] = position_list[i]
}

spe_raw$position <- factor(spe_raw$position, levels = c("TL", "TR", "BL", "BR"))

## now arrange() to get the correct order
spe_raw <- arrange(spe_raw, brnum, position)

brains = unique(spe_raw$brnum)
pdf(here("plots", "02_build_spe", "rearrangedAllSamples.pdf"), width = 8, height = 10)

for (i in seq_along(brains)){
speb <- spe_raw[, (colData(spe_raw)$in_tissue & colData(spe_raw)$brnum == brains[i])]
samples <- unique(speb$sample_id)
samples

if (length(samples) == 1){
    p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = TRUE)
} else if (length(samples) == 2){
    p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = TRUE)
    p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = TRUE)
    grid.arrange(p1, p2, nrow = 2)
} else if (length(samples) == 3){
    p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = TRUE)
    p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = TRUE)
    p3 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = TRUE)
    grid.arrange(p1, p2, p3, nrow = 2)
} else if (length(samples) == 4){
    p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = TRUE)
    p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = TRUE)
    p3 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = TRUE)
    p4 <- plotVisium(speb[, which(speb$sample_id == samples[4])], spots = TRUE)
    grid.arrange(p1, p2, p3, p4, nrow = 2)
} else if (length(samples) == 5){
    p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = TRUE)
    p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = TRUE)
    p3 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = TRUE)
    p4 <- plotVisium(speb[, which(speb$sample_id == samples[4])], spots = TRUE)
    p5 <- plotVisium(speb[, which(speb$sample_id == samples[5])], spots = TRUE)
    grid.arrange(p1, p2, p3, p4, p5, nrow = 2)}

if (length(samples) == 1){
    x = speb[, which(speb$sample_id == samples[1])]
    df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
    a1 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
} else if (length(samples) == 2){
    x = speb[, which(speb$sample_id == samples[1])]
    df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
    a1 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
    x = speb[, which(speb$sample_id == samples[2])]
    df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
    a2 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
    grid.arrange(a1, a2, nrow = 2)
} else if (length(samples) == 3){
    x = speb[, which(speb$sample_id == samples[1])]
    df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
    a1 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
    x = speb[, which(speb$sample_id == samples[2])]
    df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
    a2 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
    x = speb[, which(speb$sample_id == samples[3])]
    df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
    a3 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
    grid.arrange(a1, a2, a3, nrow = 2)
} else if (length(samples) == 4){
    x = speb[, which(speb$sample_id == samples[1])]
    df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
    a1 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
    x = speb[, which(speb$sample_id == samples[2])]
    df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
    a2 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
    x = speb[, which(speb$sample_id == samples[3])]
    df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
    a3 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
    x = speb[, which(speb$sample_id == samples[4])]
    df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
    a4 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
    grid.arrange(a1, a2, a3, a4, nrow = 2)
} else if (length(samples) == 5){
    x = speb[, which(speb$sample_id == samples[1])]
    df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
    a1 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
    x = speb[, which(speb$sample_id == samples[2])]
    df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
    a2 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
    x = speb[, which(speb$sample_id == samples[3])]
    df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
    a3 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
    x = speb[, which(speb$sample_id == samples[4])]
    df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
    a4 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
    x = speb[, which(speb$sample_id == samples[5])]
    df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
    a5 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
    grid.arrange(a1, a2, a3, a4, a5, nrow = 2)}
}

dev.off()

samples = unique(spe_raw$sample_id)
samples
# [1] "V10B01-086_D1" "V11U08-081_C1" "V10B01-086_C1" "V11U08-081_D1"
# [5] "V11L05-333_A1" "V11L05-333_B1" "V11L05-333_C1" "V11L05-333_D1"
# [9] "V10B01-085_B1" "V10B01-085_A1" "V10B01-085_D1" "V10B01-085_C1"
# [13] "V10B01-086_A1" "V10B01-086_B1" "V11L05-335_C1" "V11L05-335_B1"
# [17] "V11L05-335_A1" "V11U08-084_A1" "V11U08-084_B1" "V11U08-084_C1"
# [21] "V11U08-084_D1" "V11A20-297_C1" "V11A20-297_D1" "V11L05-335_D1"
# [25] "V11A20-297_A1" "V11A20-297_B1" "V11U08-081_A1" "V11U08-081_B1"
# [29] "V11L05-336_A1" "V11L05-336_B1" "V11L05-336_C1" "V11L05-336_D1"
# [33] "V12F14-051_C1" "V12F14-051_D1" "V12F14-051_A1" "V12F14-051_B1"

# Rotations by sample:
angle_list <- c(270, 90, 0, 0,
                0, 0, 0, 0,
                180, 90, 180, 180, 
                0, 0,
                0, 0, 0, 
                0, 0, 0, 0,
                0, 0, 0, 0, 0,
                90, 90, 
                0, 0, 0, 0,
                0,0,0,0)
angle_list_array <- c(270, 90, 0, 0,
                0, 0, 0, 0,
                180, 90, 180, 180, 
                0, 0,
                270, 180, 270, 
                0, 0, 0, 0,
                0, 0, 180, 0, 0,
                180, 180, 
                0, 0, 0, 0,
                0,0,0,90)
source(file = here::here("code", "02_build_spe", "transform_spe.R"))
source(file = here::here("code", "02_build_spe", "transform_spe_array.R"))
for (i in seq_along(angle_list)) {
    id <- samples[i]
    x <- trans_geom(spe_raw, sample_id = id, degrees = angle_list[i])
    x <- trans_spe_array(x,  degrees = angle_list_array[i])
    if (i == 1) {
        spe <- x
    } else {
        spe <- cbind(spe, x)
    }

    spe <- rotateImg(spe, sample_id = id, image_id = TRUE, angle_list[i])
   }


pdf(file = here::here("plots", "02_build_spe", "referenceMapping_allSamples.pdf"), h = 10, w = 20)
for (i in seq_along(samples)) {
    id <- samples[i]
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

brains = unique(spe$brnum)
pdf(here("plots", "02_build_spe", "rearrangedTransformedAllSamples.pdf"), width = 8, height = 10)

for (i in seq_along(brains)){
    speb <- spe[, (colData(spe)$in_tissue & colData(spe)$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    samples
    
    if (length(samples) == 1){
        p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = TRUE)
    } else if (length(samples) == 2){
        p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = TRUE)
        p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = TRUE)
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = TRUE)
        p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = TRUE)
        p3 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = TRUE)
        grid.arrange(p1, p2, p3, nrow = 2)
    } else if (length(samples) == 4){
        p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = TRUE)
        p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = TRUE)
        p3 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = TRUE)
        p4 <- plotVisium(speb[, which(speb$sample_id == samples[4])], spots = TRUE)
        grid.arrange(p1, p2, p3, p4, nrow = 2)
    } else if (length(samples) == 5){
        p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = TRUE)
        p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = TRUE)
        p3 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = TRUE)
        p4 <- plotVisium(speb[, which(speb$sample_id == samples[4])], spots = TRUE)
        p5 <- plotVisium(speb[, which(speb$sample_id == samples[5])], spots = TRUE)
        grid.arrange(p1, p2, p3, p4, p5, nrow = 2)}
    if (length(samples) == 1){
        x = speb[, which(speb$sample_id == samples[1])]
        df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
        a1 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
    } else if (length(samples) == 2){
        x = speb[, which(speb$sample_id == samples[1])]
        df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
        a1 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
        x = speb[, which(speb$sample_id == samples[2])]
        df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
        a2 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
        grid.arrange(a1, a2, nrow = 2)
    } else if (length(samples) == 3){
        x = speb[, which(speb$sample_id == samples[1])]
        df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
        a1 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
        x = speb[, which(speb$sample_id == samples[2])]
        df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
        a2 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
        x = speb[, which(speb$sample_id == samples[3])]
        df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
        a3 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
        grid.arrange(a1, a2, a3, nrow = 2)
    } else if (length(samples) == 4){
        x = speb[, which(speb$sample_id == samples[1])]
        df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
        a1 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
        x = speb[, which(speb$sample_id == samples[2])]
        df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
        a2 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
        x = speb[, which(speb$sample_id == samples[3])]
        df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
        a3 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
        x = speb[, which(speb$sample_id == samples[4])]
        df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
        a4 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
        grid.arrange(a1, a2, a3, a4, nrow = 2)
    } else if (length(samples) == 5){
        x = speb[, which(speb$sample_id == samples[1])]
        df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
        a1 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
        x = speb[, which(speb$sample_id == samples[2])]
        df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
        a2 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
        x = speb[, which(speb$sample_id == samples[3])]
        df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
        a3 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
        x = speb[, which(speb$sample_id == samples[4])]
        df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
        a4 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
        x = speb[, which(speb$sample_id == samples[5])]
        df_raw <- cbind.data.frame(colData(x), spatialCoords(x))
        a5 <- ggplot(df_raw, aes(x = array_row, y = array_col, color = sample_id)) + geom_point(size = 1) + coord_fixed() + theme_bw()
        grid.arrange(a1, a2, a3, a4, a5, nrow = 2)}
}

dev.off()


#########

## Size in Gb
lobstr::obj_size(spe)
# 6.16 GB
dim(spe)
# [1] 36601 179712
save(spe, file = here::here("processed-data", "02_build_spe", "spe_transform.Rdata"))
