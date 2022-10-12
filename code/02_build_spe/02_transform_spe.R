
setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')

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
})

load(file = here::here("processed-data", "02_build_spe", "spe_raw.Rdata"))
spe_raw = spe
dim(spe_raw)
# [1]  36601 159744

sampleID <- unique(spe_raw$sample_id)
sampleID
# [1] "V10B01-086_C1" "V10B01-086_D1" "V11U08-081_C1" "V11U08-081_D1"
# [5] "V11L05-333_A1" "V11L05-333_B1" "V11L05-333_C1" "V11L05-333_D1"
# [9] "V10B01-085_A1" "V10B01-085_B1" "V10B01-085_C1" "V10B01-085_D1"
# [13] "V10B01-086_A1" "V10B01-086_B1" "V11L05-335_A1" "V11L05-335_B1"
# [17] "V11L05-335_C1" "V11U08-084_A1" "V11U08-084_B1" "V11U08-084_C1"
# [21] "V11U08-084_D1" "V11A20-297_A1" "V11A20-297_B1" "V11A20-297_C1"
# [25] "V11A20-297_D1" "V11L05-335_D1" "V11U08-081_A1" "V11U08-081_B1"
# [29] "V11L05-336_A1" "V11L05-336_B1" "V11L05-336_C1" "V11L05-336_D1"

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
pdf(here("plots", "02_build_spe", "prearrangedSamples.pdf"), width = 8, height = 10)

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

#                   Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325 Br8492 Br8667
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

samples = unique(spe_raw$sample_id)
# Br2743 [1] "V10B01-086_C1" "V11U08-081_C1" "V10B01-086_D1" "V11U08-081_D1"
# Br3942 [5] "V11L05-333_A1" "V11L05-333_B1" "V11L05-333_C1" "V11L05-333_D1"
# Br6423 [9] "V10B01-085_A1" "V10B01-085_B1" "V10B01-085_C1" "V10B01-085_D1"
# Br6432 [13] "V10B01-086_A1" "V10B01-086_B1" 
# Br6471 [15] "V11L05-335_A1" "V11L05-335_B1" "V11L05-335_C1" 
# Br6522 [18] "V11U08-084_A1" "V11U08-084_B1" "V11U08-084_C1" "V11U08-084_D1" 
# Br8325 [22] "V11A20-297_A1" "V11A20-297_B1" "V11A20-297_C1" "V11A20-297_D1"  "V11L05-335_D1" 
# Br8492 [27] "V11U08-081_A1" "V11U08-081_B1"
# Br8667 [29] "V11L05-336_A1" "V11L05-336_B1" "V11L05-336_C1" "V11L05-336_D1"

position_list <- c("TR", "BL", "TL", "BR",
                   "TL", "TR", "BL", "BR",
                   "TR", "TL", "BR", "BL",
                   "TL", "BL",
                   "BL", "TR", "TL",
                   "TL", "TR", "BL", "BR",
                   "BL", "BR", "TL", "TR", "TR",
                   "TL", "BL",
                   "TL", "TR", "BL", "BR")

for (i in seq_along(samples)){
    spe_raw$position[colData(spe_raw)$sample_id == samples[i]] = position_list[i]
}

spe_raw$position <- factor(spe_raw$position, levels = c("TL", "TR", "BL", "BR"))

## now arrange() to get the correct order
spe_raw <- arrange(spe_raw, brnum, position)

brains = unique(spe_raw$brnum)
pdf(here("plots", "02_build_spe", "rearrangedSamples.pdf"), width = 8, height = 10)

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
# [1] "V10B01-086_D1" "V10B01-086_C1" "V11U08-081_C1" "V11U08-081_D1"
# [5] "V11L05-333_A1" "V11L05-333_B1" "V11L05-333_C1" "V11L05-333_D1"
# [9] "V10B01-085_B1" "V10B01-085_A1" "V10B01-085_D1" "V10B01-085_C1"
# [13] "V10B01-086_A1" "V10B01-086_B1" "V11L05-335_C1" "V11L05-335_B1"
# [17] "V11L05-335_A1" "V11U08-084_A1" "V11U08-084_B1" "V11U08-084_C1"
# [21] "V11U08-084_D1" "V11A20-297_C1" "V11A20-297_D1" "V11L05-335_D1"
# [25] "V11A20-297_A1" "V11A20-297_B1" "V11U08-081_A1" "V11U08-081_B1"
# [29] "V11L05-336_A1" "V11L05-336_B1" "V11L05-336_C1" "V11L05-336_D1"

# Rotations by sample:
angle_list <- c(270, 90, 0, 0,
                0, 0, 0, 0,
                180, 90, 180, 180, 
                0, 0,
                0, 0, 0, 
                0, 0, 0, 0,
                0, 0, 0, 0, 0,
                90, 90, 
                0, 0, 0, 0)
angle_list_array <- c(270, 90, 0, 0,
                0, 0, 0, 0,
                180, 90, 180, 180, 
                0, 0,
                270, 180, 270, 
                0, 0, 0, 0,
                0, 0, 180, 0, 0,
                180, 180, 
                0, 0, 0, 0)
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


pdf(file = here::here("plots", "02_build_spe", "referenceMapping.pdf"), h = 10, w = 20)
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
pdf(here("plots", "02_build_spe", "rearrangedTransformedSamples.pdf"), width = 8, height = 10)

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
# 5.53 GB
dim(spe)
# [1] 36601 159744
save(spe, file = here::here("processed-data", "02_build_spe", "spe_transform.Rdata"))

## Reproducibility information
# Sys.time
# "2022-10-05 12:48:52 EDT"
# proc.time()
# user   system  elapsed 
# 786.370  131.882 6005.361 
# 
# options(width = 120)
# session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R Under development (unstable) (2021-11-06 r81149)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-10-05
# pandoc   2.11.0.4 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                  * version   date (UTC) lib source
# AnnotationDbi              1.58.0    2022-04-26 [1] Bioconductor
# AnnotationHub              3.4.0     2022-04-26 [1] Bioconductor
# assertthat                 0.2.1     2019-03-21 [2] CRAN (R 4.1.0)
# attempt                    0.3.1     2020-05-03 [1] CRAN (R 4.2.0)
# beachmat                   2.12.0    2022-04-26 [1] Bioconductor
# beeswarm                   0.4.0     2021-06-01 [2] CRAN (R 4.2.0)
# benchmarkme                1.0.7     2021-03-21 [1] CRAN (R 4.2.0)
# benchmarkmeData            1.0.4     2020-04-23 [1] CRAN (R 4.2.0)
# Biobase                  * 2.56.0    2022-04-26 [1] Bioconductor
# BiocFileCache              2.4.0     2022-04-26 [1] Bioconductor
# BiocGenerics             * 0.42.0    2022-04-26 [1] Bioconductor
# BiocIO                     1.6.0     2022-04-26 [1] Bioconductor
# BiocManager                1.30.18   2022-05-18 [2] CRAN (R 4.2.0)
# BiocNeighbors              1.14.0    2022-04-26 [1] Bioconductor
# BiocParallel               1.30.3    2022-06-05 [1] Bioconductor
# BiocSingular               1.12.0    2022-04-26 [1] Bioconductor
# BiocVersion                3.15.2    2022-03-29 [2] Bioconductor
# Biostrings                 2.64.0    2022-04-26 [1] Bioconductor
# bit                        4.0.4     2020-08-04 [2] CRAN (R 4.1.0)
# bit64                      4.0.5     2020-08-30 [2] CRAN (R 4.1.0)
# bitops                     1.0-7     2021-04-24 [2] CRAN (R 4.2.0)
# blob                       1.2.3     2022-04-10 [2] CRAN (R 4.2.0)
# bslib                      0.4.0     2022-07-16 [2] CRAN (R 4.2.0)
# cachem                     1.0.6     2021-08-19 [2] CRAN (R 4.2.0)
# cli                        3.4.1     2022-09-23 [2] CRAN (R 4.2.0)
# codetools                  0.2-18    2020-11-04 [3] CRAN (R 4.2.0)
# colorout                 * 1.2-2     2022-03-10 [1] Github (jalvesaq/colorout@79931fd)
# colorspace                 2.0-3     2022-02-21 [2] CRAN (R 4.2.0)
# config                     0.3.1     2020-12-17 [1] CRAN (R 4.2.0)
# cowplot                    1.1.1     2020-12-30 [2] CRAN (R 4.2.0)
# crayon                     1.5.2     2022-09-29 [2] CRAN (R 4.2.0)
# curl                       4.3.2     2021-06-23 [2] CRAN (R 4.2.0)
# data.table                 1.14.2    2021-09-27 [2] CRAN (R 4.2.0)
# DBI                        1.1.3     2022-06-18 [2] CRAN (R 4.2.0)
# dbplyr                     2.2.1     2022-06-27 [2] CRAN (R 4.2.0)
# DelayedArray               0.22.0    2022-04-26 [1] Bioconductor
# DelayedMatrixStats         1.18.0    2022-04-26 [1] Bioconductor
# desc                       1.4.2     2022-09-08 [2] CRAN (R 4.2.0)
# digest                     0.6.29    2021-12-01 [2] CRAN (R 4.2.0)
# doParallel                 1.0.17    2022-02-07 [2] CRAN (R 4.2.0)
# dotCall64                  1.0-2     2022-10-03 [2] CRAN (R 4.2.0)
# dplyr                      1.0.10    2022-09-01 [2] CRAN (R 4.2.0)
# dqrng                      0.3.0     2021-05-01 [2] CRAN (R 4.2.0)
# DropletUtils               1.16.0    2022-04-26 [2] Bioconductor
# DT                         0.25      2022-09-12 [2] CRAN (R 4.2.0)
# edgeR                      3.38.1    2022-05-15 [1] Bioconductor
# ellipsis                   0.3.2     2021-04-29 [2] CRAN (R 4.2.0)
# ExperimentHub              2.4.0     2022-04-26 [1] Bioconductor
# fansi                      1.0.3     2022-03-24 [2] CRAN (R 4.2.0)
# farver                     2.1.1     2022-07-06 [2] CRAN (R 4.2.0)
# fastmap                    1.1.0     2021-01-25 [2] CRAN (R 4.1.0)
# fields                     14.1      2022-08-12 [2] CRAN (R 4.2.0)
# filelock                   1.0.2     2018-10-05 [2] CRAN (R 4.1.0)
# foreach                    1.5.2     2022-02-02 [2] CRAN (R 4.2.0)
# fs                         1.5.2     2021-12-08 [2] CRAN (R 4.2.0)
# generics                   0.1.3     2022-07-05 [2] CRAN (R 4.2.0)
# GenomeInfoDb             * 1.32.2    2022-05-15 [1] Bioconductor
# GenomeInfoDbData           1.2.8     2022-04-16 [2] Bioconductor
# GenomicAlignments          1.32.0    2022-04-26 [1] Bioconductor
# GenomicRanges            * 1.48.0    2022-04-26 [1] Bioconductor
# ggbeeswarm                 0.6.0     2017-08-07 [2] CRAN (R 4.2.0)
# ggplot2                  * 3.3.6     2022-05-03 [2] CRAN (R 4.2.0)
# ggrepel                    0.9.1     2021-01-15 [2] CRAN (R 4.1.0)
# ggside                     0.2.1     2022-07-20 [2] CRAN (R 4.2.0)
# ggspavis                 * 1.2.0     2022-04-26 [2] Bioconductor
# glue                       1.6.2     2022-02-24 [2] CRAN (R 4.2.0)
# golem                      0.3.2     2022-03-04 [1] CRAN (R 4.2.0)
# gridExtra                * 2.3       2017-09-09 [2] CRAN (R 4.1.0)
# gtable                     0.3.1     2022-09-01 [2] CRAN (R 4.2.0)
# HDF5Array                  1.24.1    2022-06-02 [1] Bioconductor
# here                     * 1.0.1     2020-12-13 [1] CRAN (R 4.2.0)
# htmltools                  0.5.3     2022-07-18 [2] CRAN (R 4.2.0)
# htmlwidgets                1.5.4     2021-09-08 [2] CRAN (R 4.2.0)
# httpuv                     1.6.6     2022-09-08 [2] CRAN (R 4.2.0)
# httr                       1.4.4     2022-08-17 [2] CRAN (R 4.2.0)
# interactiveDisplayBase     1.34.0    2022-04-26 [1] Bioconductor
# IRanges                  * 2.30.0    2022-04-26 [1] Bioconductor
# irlba                      2.3.5.1   2022-10-03 [2] CRAN (R 4.2.0)
# iterators                  1.0.14    2022-02-05 [2] CRAN (R 4.2.0)
# jquerylib                  0.1.4     2021-04-26 [2] CRAN (R 4.2.0)
# jsonlite                   1.8.2     2022-10-02 [2] CRAN (R 4.2.0)
# KEGGREST                   1.36.0    2022-04-26 [1] Bioconductor
# knitr                      1.40      2022-08-24 [2] CRAN (R 4.2.0)
# labeling                   0.4.2     2020-10-20 [2] CRAN (R 4.1.0)
# later                      1.3.0     2021-08-18 [2] CRAN (R 4.2.0)
# lattice                    0.20-45   2021-09-22 [3] CRAN (R 4.2.0)
# lazyeval                   0.2.2     2019-03-15 [2] CRAN (R 4.1.0)
# lifecycle                  1.0.2     2022-09-09 [2] CRAN (R 4.2.0)
# limma                      3.52.2    2022-06-19 [1] Bioconductor
# lobstr                   * 1.1.2     2022-06-22 [2] CRAN (R 4.2.0)
# locfit                     1.5-9.6   2022-07-11 [2] CRAN (R 4.2.0)
# magick                     2.7.3     2021-08-18 [2] CRAN (R 4.2.0)
# magrittr                   2.0.3     2022-03-30 [2] CRAN (R 4.2.0)
# maps                       3.4.0     2021-09-25 [2] CRAN (R 4.2.0)
# Matrix                     1.5-1     2022-09-13 [2] CRAN (R 4.2.0)
# MatrixGenerics           * 1.8.1     2022-06-26 [1] Bioconductor
# matrixStats              * 0.62.0    2022-04-19 [2] CRAN (R 4.2.0)
# memoise                    2.0.1     2021-11-26 [2] CRAN (R 4.2.0)
# mime                       0.12      2021-09-28 [2] CRAN (R 4.2.0)
# munsell                    0.5.0     2018-06-12 [2] CRAN (R 4.1.0)
# paletteer                  1.4.0     2021-07-20 [1] CRAN (R 4.2.0)
# pillar                     1.8.1     2022-08-19 [2] CRAN (R 4.2.0)
# pkgconfig                  2.0.3     2019-09-22 [2] CRAN (R 4.1.0)
# pkgload                    1.3.0     2022-06-27 [2] CRAN (R 4.2.0)
# plotly                     4.10.0    2021-10-09 [2] CRAN (R 4.2.0)
# png                        0.1-7     2013-12-03 [2] CRAN (R 4.1.0)
# prettyunits                1.1.1     2020-01-24 [2] CRAN (R 4.1.0)
# promises                   1.2.0.1   2021-02-11 [2] CRAN (R 4.1.0)
# purrr                      0.3.4     2020-04-17 [2] CRAN (R 4.1.0)
# R.methodsS3                1.8.2     2022-06-13 [2] CRAN (R 4.2.0)
# R.oo                       1.25.0    2022-06-12 [2] CRAN (R 4.2.0)
# R.utils                    2.12.0    2022-06-28 [2] CRAN (R 4.2.0)
# R6                         2.5.1     2021-08-19 [2] CRAN (R 4.2.0)
# rappdirs                   0.3.3     2021-01-31 [2] CRAN (R 4.1.0)
# RColorBrewer               1.1-3     2022-04-03 [2] CRAN (R 4.2.0)
# Rcpp                       1.0.9     2022-07-08 [2] CRAN (R 4.2.0)
# RCurl                      1.98-1.9  2022-10-03 [2] CRAN (R 4.2.0)
# rematch2                   2.1.2     2020-05-01 [2] CRAN (R 4.1.0)
# restfulr                   0.0.15    2022-06-16 [2] CRAN (R 4.2.0)
# rhdf5                      2.40.0    2022-04-26 [1] Bioconductor
# rhdf5filters               1.8.0     2022-04-26 [1] Bioconductor
# Rhdf5lib                   1.18.2    2022-05-15 [1] Bioconductor
# rjson                      0.2.21    2022-01-09 [2] CRAN (R 4.2.0)
# rlang                      1.0.5     2022-08-31 [1] CRAN (R 4.2.0)
# roxygen2                   7.2.1     2022-07-18 [2] CRAN (R 4.2.0)
# rprojroot                  2.0.3     2022-04-02 [2] CRAN (R 4.2.0)
# Rsamtools                  2.12.0    2022-04-26 [1] Bioconductor
# RSQLite                    2.2.18    2022-10-04 [2] CRAN (R 4.2.0)
# rstudioapi                 0.14      2022-08-22 [2] CRAN (R 4.2.0)
# rsvd                       1.0.5     2021-04-16 [2] CRAN (R 4.2.0)
# rtracklayer              * 1.56.0    2022-04-26 [1] Bioconductor
# S4Vectors                * 0.34.0    2022-04-26 [1] Bioconductor
# sass                       0.4.2     2022-07-16 [2] CRAN (R 4.2.0)
# ScaledMatrix               1.4.0     2022-04-26 [1] Bioconductor
# scales                     1.2.1     2022-08-20 [2] CRAN (R 4.2.0)
# scater                     1.24.0    2022-04-26 [1] Bioconductor
# scuttle                    1.6.2     2022-05-15 [1] Bioconductor
# sessioninfo              * 1.2.2     2021-12-06 [2] CRAN (R 4.2.0)
# shiny                      1.7.2     2022-07-19 [2] CRAN (R 4.2.0)
# shinyWidgets               0.7.3     2022-08-31 [2] CRAN (R 4.2.0)
# SingleCellExperiment     * 1.18.0    2022-04-26 [1] Bioconductor
# spam                       2.9-1     2022-08-07 [2] CRAN (R 4.2.0)
# sparseMatrixStats          1.8.0     2022-04-26 [1] Bioconductor
# SpatialExperiment        * 1.6.0     2022-04-26 [1] Bioconductor
# spatialLIBD              * 1.9.2     2022-05-13 [1] Github (LieberInstitute/spatialLIBD@35ccde7)
# stringi                    1.7.8     2022-07-11 [2] CRAN (R 4.2.0)
# stringr                    1.4.1     2022-08-20 [2] CRAN (R 4.2.0)
# SummarizedExperiment     * 1.26.1    2022-04-29 [1] Bioconductor
# tibble                     3.1.8     2022-07-22 [2] CRAN (R 4.2.0)
# tidyr                      1.2.1     2022-09-08 [2] CRAN (R 4.2.0)
# tidyselect                 1.1.2     2022-02-21 [2] CRAN (R 4.2.0)
# tidySingleCellExperiment * 1.6.3     2022-05-22 [1] Bioconductor
# ttservice                * 0.2.2     2022-06-24 [1] CRAN (R 4.2.0)
# usethis                    2.1.6     2022-05-25 [2] CRAN (R 4.2.0)
# utf8                       1.2.2     2021-07-24 [2] CRAN (R 4.2.0)
# vctrs                      0.4.2     2022-09-29 [2] CRAN (R 4.2.0)
# vipor                      0.4.5     2017-03-22 [2] CRAN (R 4.2.0)
# viridis                    0.6.2     2021-10-13 [2] CRAN (R 4.2.0)
# viridisLite                0.4.1     2022-08-22 [2] CRAN (R 4.2.0)
# withr                      2.5.0     2022-03-03 [2] CRAN (R 4.2.0)
# xfun                       0.33      2022-09-12 [2] CRAN (R 4.2.0)
# XML                        3.99-0.11 2022-10-03 [2] CRAN (R 4.2.0)
# xml2                       1.3.3     2021-11-30 [2] CRAN (R 4.2.0)
# xtable                     1.8-4     2019-04-21 [2] CRAN (R 4.1.0)
# XVector                    0.36.0    2022-04-26 [1] Bioconductor
# yaml                       2.3.5     2022-02-21 [2] CRAN (R 4.2.0)
# zlibbioc                   1.42.0    2022-04-26 [1] Bioconductor
# 
# [1] /users/mtippani/R/devel
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/R/devel/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/R/devel/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
