
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

# Rotations by sample:

# V10B01-085_A1= 90, TR
# V10B01-085_B1= 180, TL
# V10B01-085_C1= 180, BR
# V10B01-085_D1= 180, BL
# V10B01-086_A1= 0, TL
# V10B01-086_B1= 0, BL
# V10B01-086_C1= 0, TR
# V10B01-086_D1= 180, BR
# V11L05-333_A1= 0, TL
# V11L05-333_B1= 0, TR
# V11L05-333_C1= 0, BL
# V11L05-333_D1= 0, BR
# V11L05-335_A1= 0, BL
# V11L05-335_B1= 0, TR
# V11L05-335_C1= 0, TL
# V11L05-335_D1= 0, BR
# V11U08-084_A1= 0, TL
# V11U08-084_B1= 0, TR
# V11U08-084_C1= 0, BL
# V11U08-084_D1= 0, BR
# V11A20-297_A1= 0, BL
# V11A20-297_B1= 0, BR
# V11A20-297_C1= 0, TL
# V11A20-297_D1= 0, TR
# V11U08-081_A1= 0, TL
# V11U08-081_B1= 0, TR
# V11U08-081_C1= 0, BL
# V11U08-081_D1= 0, BR
# V11L05-336_A1= 0, BR
# V11L05-336_B1= 0, BL
# V11L05-336_C1= 0, TR
# V11L05-336_D1= 0, TL

sampleID <- unique(spe_raw$sample_id)
# [1] "V10B01-086_C1" "V10B01-086_D1" "V11U08-081_C1" "V11U08-081_D1"
# [5] "V11L05-333_A1" "V11L05-333_B1" "V11L05-333_C1" "V11L05-333_D1"
# [9] "V10B01-085_A1" "V10B01-085_B1" "V10B01-085_C1" "V10B01-085_D1"
# [13] "V10B01-086_A1" "V10B01-086_B1" "V11L05-335_A1" "V11L05-335_B1"
# [17] "V11L05-335_C1" "V11U08-084_A1" "V11U08-084_B1" "V11U08-084_C1"
# [21] "V11U08-084_D1" "V11A20-297_A1" "V11A20-297_B1" "V11A20-297_C1"
# [25] "V11A20-297_D1" "V11L05-335_D1" "V11U08-081_A1" "V11U08-081_B1"
# [29] "V11L05-336_A1" "V11L05-336_B1" "V11L05-336_C1" "V11L05-336_D1"

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

unique(spe$sample_id[which(spe$position == "N")])
spe$position[which(spe$sample_id == "V11L05-335_D1")] <- "BR"

unique(spe$sample_id[which(spe$position == "")])
spe$position[which(spe$sample_id == "V10B01-085_A1")] <- "TR"
spe$position[which(spe$sample_id == "V10B01-085_B1")] <- "TL"
spe$position[which(spe$sample_id == "V10B01-085_C1")] <- "BR"
spe$position[which(spe$sample_id == "V10B01-085_D1")] <- "BL"
spe$position <- factor(spe$position, levels = c("TL", "TR", "BR", "BL"))

unique(spe$position)

slide_order <- unique(spe$slide)
sample_order <- unlist(sapply(slide_order, function(i) {
    (unique(spe$sample_id)[grepl(i, unique(spe$sample_id))])
}))
sample_order

# [,1]            [,2]            [,3]            [,4]
# [1,] "V10B01-086_C1" "V11U08-081_C1" "V11L05-333_A1" "V10B01-085_A1"
# [2,] "V10B01-086_D1" "V11U08-081_D1" "V11L05-333_B1" "V10B01-085_B1"
# [3,] "V10B01-086_A1" "V11U08-081_A1" "V11L05-333_C1" "V10B01-085_C1"
# [4,] "V10B01-086_B1" "V11U08-081_B1" "V11L05-333_D1" "V10B01-085_D1"
# [,5]            [,6]            [,7]            [,8]
# [1,] "V11L05-335_A1" "V11U08-084_A1" "V11A20-297_A1" "V11L05-336_A1"
# [2,] "V11L05-335_B1" "V11U08-084_B1" "V11A20-297_B1" "V11L05-336_B1"
# [3,] "V11L05-335_C1" "V11U08-084_C1" "V11A20-297_C1" "V11L05-336_C1"
# [4,] "V11L05-335_D1" "V11U08-084_D1" "V11A20-297_D1" "V11L05-336_D1"

## code for rearranging spe for better Shiny app visualization
library(tidySingleCellExperiment)
## change colnames(spe) and rownames(spatialCoords(spe)) to spe$key
colnames(spe) <- spe$key
rownames(spatialCoords(spe)) <- spe$key
## now arrange() to get the correct order
spe <- arrange(spe, slide, array)


Newsample_order <- unlist(sapply(slide_order, function(i) {
    (unique(temp$sample_id)[grepl(i, unique(temp$sample_id))])
}))
Newsample_order

# [,1]            [,2]            [,3]            [,4]
# [1,] "V10B01-086_A1" "V11U08-081_C1" "V11L05-333_A1" "V10B01-085_B1"
# [2,] "V10B01-086_C1" "V11U08-081_A1" "V11L05-333_B1" "V10B01-085_A1"
# [3,] "V10B01-086_D1" "V11U08-081_D1" "V11L05-333_D1" "V10B01-085_C1"
# [4,] "V10B01-086_B1" "V11U08-081_B1" "V11L05-333_C1" "V10B01-085_D1"
# [,5]            [,6]            [,7]            [,8]
# [1,] "V11L05-335_C1" "V11U08-084_A1" "V11A20-297_C1" "V11L05-336_D1"
# [2,] "V11L05-335_B1" "V11U08-084_B1" "V11A20-297_D1" "V11L05-336_C1"
# [3,] "V11L05-335_A1" "V11U08-084_D1" "V11A20-297_B1" "V11L05-336_A1"
# [4,] "V11L05-335_D1" "V11U08-084_C1" "V11A20-297_A1" "V11L05-336_B1"

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)
length(no_expr)
# [1] 6242
length(no_expr) / nrow(spe) * 100
# [1] 17.05418
spe <- spe[-no_expr, ]

## Now drop the spots outside the tissue
spe <- spe[, colData(spe)$in_tissue]

## Size in Gb
lobstr::obj_size(spe)
# 4.856146 B
dim(spe)
# [1] 30359 137446

## Remove spots without counts
if (any(colSums(counts(spe)) == 0)) {
    message("removing spots without counts for spe")
    spe <- spe[, -which(colSums(counts(spe)) == 0)]
    dim(spe)
}

lobstr::obj_size(spe)
# 4.856144 B
dim(spe)
# [1]  30359 137440

save(spe, file = here::here("processed-data", "02_build_spe", "spe_basic.Rdata"))
