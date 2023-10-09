setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
  library("here")
  library("sessioninfo")
  library("SpatialExperiment")
  library("spatialLIBD")
  library("Polychrome")
  library("ggplot2")
})

offset_check <- function(spe){
speo = spe

# The spatial locations for each sample need to be offset so that spots of different samples are not neighbors.
# summary(colData(spe)$array_row)
# summary(colData(spe)$array_col)

brains = unique(speo$brnum)
brains[1]
#Br2743

x = speo[,colData(speo)$brnum == brains[1]]
samples = unique(x$sample_id)
samples
# [1] "V10B01-086_D1" "V10B01-086_C1" "V11U08-081_C1" "V11U08-081_D1"

row <- x$array_row
col <- x$array_col

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 128
ix = colData(x)$sample_id == samples[2]
col[ix] <- col[ix] - 4 

ix = colData(x)$sample_id == samples[3]
col[ix] <- col[ix] + 100 
ix = colData(x)$sample_id == samples[4]
row[ix] <- row[ix] + 77 
ix = colData(x)$sample_id == samples[4]
col[ix] <- col[ix] + 101 

colData(x)$row <- row
colData(x)$col <- col

spe = x

brains[2]
#Br3942

x = speo[,colData(speo)$brnum == brains[2]]
samples = unique(x$sample_id)
samples
# [1] "V11L05-333_A1" "V11L05-333_B1" "V11L05-333_C1" "V11L05-333_D1"

row <- x$array_row
col <- x$array_col

ix = colData(x)$sample_id == samples[1]
col[ix] <- col[ix] + 398

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 78
ix = colData(x)$sample_id == samples[2]
col[ix] <- col[ix] + 398

ix = colData(x)$sample_id == samples[3]
col[ix] <- col[ix] + 270 

ix = colData(x)$sample_id == samples[4]
row[ix] <- row[ix] + 78 
ix = colData(x)$sample_id == samples[4]
col[ix] <- col[ix] + 270 

colData(x)$row <- row
colData(x)$col <- col

spe = cbind(spe,x)

brains[3]
#Br6423

x = speo[,colData(speo)$brnum == brains[3]]
samples = unique(x$sample_id)
samples
# [1] "V10B01-085_B1" "V10B01-085_A1" "V10B01-085_D1" "V10B01-085_C1"

row <- x$array_row
col <- x$array_col

ix = colData(x)$sample_id == samples[1]
row[ix] <- row[ix] + 13
ix = colData(x)$sample_id == samples[1]
col[ix] <- col[ix] + 659

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 76
ix = colData(x)$sample_id == samples[2]
col[ix] <- col[ix] + 666

ix = colData(x)$sample_id == samples[3]
col[ix] <- col[ix] + 570 

ix = colData(x)$sample_id == samples[4]
row[ix] <- row[ix] + 75
ix = colData(x)$sample_id == samples[4]
col[ix] <- col[ix] + 539

colData(x)$row <- row
colData(x)$col <- col

spe = cbind(spe,x)

brains[4]
#Br6432

x = speo[,colData(speo)$brnum == brains[4]]
samples = unique(x$sample_id)
samples
# [1] "V10B01-086_A1" "V10B01-086_B1"

row <- x$array_row
col <- x$array_col

ix = colData(x)$sample_id == samples[1]
row[ix] <- row[ix] + 300
ix = colData(x)$sample_id == samples[1]
col[ix] <- col[ix] + 128

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 300


colData(x)$row <- row
colData(x)$col <- col

spe = cbind(spe,x)

brains[5]
#Br6471

x = speo[,colData(speo)$brnum == brains[5]]
samples = unique(x$sample_id)
samples
# [1] "V11L05-335_C1" "V11L05-335_B1" "V11L05-335_A1"

row <- x$array_row
col <- x$array_col

ix = colData(x)$sample_id == samples[1]
row[ix] <- row[ix] + 300
ix = colData(x)$sample_id == samples[1]
col[ix] <- col[ix] + 380

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 428
ix = colData(x)$sample_id == samples[2]
col[ix] <- col[ix] + 399

ix = colData(x)$sample_id == samples[3]
row[ix] <- row[ix] + 300
ix = colData(x)$sample_id == samples[3]
col[ix] <- col[ix] + 300

colData(x)$row <- row
colData(x)$col <- col

spe = cbind(spe,x)

brains[6]
#Br6522

x = speo[,colData(speo)$brnum == brains[6]]
samples = unique(x$sample_id)
samples
# [1] "V11U08-084_A1" "V11U08-084_B1" "V11U08-084_C1" "V11U08-084_D1"

row <- x$array_row
col <- x$array_col

ix = colData(x)$sample_id == samples[1]
row[ix] <- row[ix] + 295
ix = colData(x)$sample_id == samples[1]
col[ix] <- col[ix] + 731

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 374
ix = colData(x)$sample_id == samples[2]
col[ix] <- col[ix] + 700

ix = colData(x)$sample_id == samples[3]
row[ix] <- row[ix] + 295
ix = colData(x)$sample_id == samples[3]
col[ix] <- col[ix] + 602

ix = colData(x)$sample_id == samples[4]
row[ix] <- row[ix] + 374
ix = colData(x)$sample_id == samples[4]
col[ix] <- col[ix] + 571

colData(x)$row <- row
colData(x)$col <- col

spe = cbind(spe,x)

brains[7]
#Br8325

x = speo[,colData(speo)$brnum == brains[7]]
samples = unique(x$sample_id)
samples
# [1] "V11U08-084_A1" "V11U08-084_B1" "V11U08-084_C1" "V11U08-084_D1"

row <- x$array_row
col <- x$array_col

ix = colData(x)$sample_id == samples[1]
row[ix] <- row[ix] + 477
ix = colData(x)$sample_id == samples[1]
col[ix] <- col[ix] + 125

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 555
ix = colData(x)$sample_id == samples[2]
col[ix] <- col[ix] + 135

ix = colData(x)$sample_id == samples[3]
row[ix] <- row[ix] + 635
ix = colData(x)$sample_id == samples[3]
col[ix] <- col[ix] + 135

ix = colData(x)$sample_id == samples[4]
row[ix] <- row[ix] + 475

ix = colData(x)$sample_id == samples[5]
row[ix] <- row[ix] + 555
ix = colData(x)$sample_id == samples[5]
col[ix] <- col[ix] + 7

colData(x)$row <- row
colData(x)$col <- col

spe = cbind(spe,x)

brains[8]
#Br8492

x = speo[,colData(speo)$brnum == brains[8]]
samples = unique(x$sample_id)
samples
# [1] "V11U08-084_A1" "V11U08-084_B1" "V11U08-084_C1" "V11U08-084_D1"

row <- x$array_row
col <- x$array_col

ix = colData(x)$sample_id == samples[1]
row[ix] <- row[ix] + 550
ix = colData(x)$sample_id == samples[1]
col[ix] <- col[ix] + 412

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 565
ix = colData(x)$sample_id == samples[2]
col[ix] <- col[ix] + 284

colData(x)$row <- row
colData(x)$col <- col

spe = cbind(spe,x)

brains[9]
#Br8667

x = speo[,colData(speo)$brnum == brains[9]]
samples = unique(x$sample_id)
samples
# [1] "V11U08-084_A1" "V11U08-084_B1" "V11U08-084_C1" "V11U08-084_D1"

row <- x$array_row
col <- x$array_col

ix = colData(x)$sample_id == samples[1]
row[ix] <- row[ix] + 550
ix = colData(x)$sample_id == samples[1]
col[ix] <- col[ix] + 740

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 628
ix = colData(x)$sample_id == samples[2]
col[ix] <- col[ix] + 740

ix = colData(x)$sample_id == samples[3]
row[ix] <- row[ix] + 550
ix = colData(x)$sample_id == samples[3]
col[ix] <- col[ix] + 595

ix = colData(x)$sample_id == samples[4]
row[ix] <- row[ix] + 626
ix = colData(x)$sample_id == samples[4]
col[ix] <- col[ix] + 580

colData(x)$row <- row
colData(x)$col <- col

spe = cbind(spe,x)

brains[9]
#Br8667

x = speo[,colData(speo)$brnum == brains[10]]
samples = unique(x$sample_id)
samples
# [1] "V11U08-084_A1" "V11U08-084_B1" "V11U08-084_C1" "V11U08-084_D1"

row <- x$array_row
col <- x$array_col

ix = colData(x)$sample_id == samples[1]
row[ix] <- row[ix] + 750
ix = colData(x)$sample_id == samples[1]
col[ix] <- col[ix] + 940

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 828
ix = colData(x)$sample_id == samples[2]
col[ix] <- col[ix] + 910

ix = colData(x)$sample_id == samples[3]
row[ix] <- row[ix] + 750
ix = colData(x)$sample_id == samples[3]
col[ix] <- col[ix] + 810

ix = colData(x)$sample_id == samples[4]
row[ix] <- row[ix] + 826
ix = colData(x)$sample_id == samples[4]
col[ix] <- col[ix] + 850

colData(x)$row <- row
colData(x)$col <- col

spe = cbind(spe,x)

#vspg Br3942
x = speo[,colData(speo)$brnum == brains[11]]
samples = unique(x$sample_id)
samples


row <- x$array_row
col <- x$array_col

ix = colData(x)$sample_id == samples[1]
col[ix] <- col[ix] + 1007

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 85
ix = colData(x)$sample_id == samples[2]
col[ix] <- col[ix] + 1001

ix = colData(x)$sample_id == samples[3]
col[ix] <- col[ix] + 887 

ix = colData(x)$sample_id == samples[4]
row[ix] <- row[ix] + 83
ix = colData(x)$sample_id == samples[4]
col[ix] <- col[ix] + 870 

colData(x)$row <- row
colData(x)$col <- col

spe = cbind(spe,x)

#vspg Br8325
x = speo[,colData(speo)$brnum == brains[12]]
samples = unique(x$sample_id)
samples


row <- x$array_row
col <- x$array_col

ix = colData(x)$sample_id == samples[1]
row[ix] <- row[ix] + 250
ix = colData(x)$sample_id == samples[1]
col[ix] <- col[ix] + 1055

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 338
ix = colData(x)$sample_id == samples[2]
col[ix] <- col[ix] + 1048

ix = colData(x)$sample_id == samples[3]
row[ix] <- row[ix] + 252
ix = colData(x)$sample_id == samples[3]
col[ix] <- col[ix] + 920 

ix = colData(x)$sample_id == samples[4]
row[ix] <- row[ix] + 330
ix = colData(x)$sample_id == samples[4]
col[ix] <- col[ix] + 920 

colData(x)$row <- row
colData(x)$col <- col

spe = cbind(spe,x)


# df <- cbind.data.frame(colData(x), spatialCoords(x))
# ggplot(df, aes(x = row, y = col, color = sample_id)) + 
#   geom_point(size = 1) + 
#   #scale_color_manual(values = cols) + 
#   coord_fixed() + 
#   guides(color = guide_legend(override.aes = list(size = 3))) + 
#   theme_bw()

# pdf(here("plots", "06_clustering", "arrayCoord_offset_updated.pdf"), width = 8, height = 10)
# df <- cbind.data.frame(colData(spe), spatialCoords(spe))
# ggplot(df, aes(x = row, y = col, color = sample_id)) +
#   geom_point(size = 1) +
#   #scale_color_manual(values = cols) +
#   coord_fixed() +
#   guides(color = guide_legend(override.aes = list(size = 3))) +
#   theme_bw()
# dev.off()

return(spe)
}
