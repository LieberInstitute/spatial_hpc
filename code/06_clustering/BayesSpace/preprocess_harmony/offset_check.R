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
brains = unique(spe$brnum)

# The spatial locations for each sample need to be offset so that spots of different samples are not neighbors.
# summary(colData(spe)$array_row)
# summary(colData(spe)$array_col)
speo = spe

brains[1]
#Br2743

x = speo[,colData(speo)$brnum == brains[1]]
samples = unique(x$sample_id)
samples
# [1] "V10B01-086_D1" "V10B01-086_C1" "V11U08-081_C1" "V11U08-081_D1"

row <- x$array_row
col <- x$array_col

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 130
ix = colData(x)$sample_id == samples[2]
col[ix] <- col[ix] - 3 

ix = colData(x)$sample_id == samples[3]
col[ix] <- col[ix] + 100 
ix = colData(x)$sample_id == samples[4]
row[ix] <- row[ix] + 78 
ix = colData(x)$sample_id == samples[4]
col[ix] <- col[ix] + 100 

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
col[ix] <- col[ix] + 400

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 78
ix = colData(x)$sample_id == samples[2]
col[ix] <- col[ix] + 400

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
col[ix] <- col[ix] + 661

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 67
ix = colData(x)$sample_id == samples[2]
col[ix] <- col[ix] + 670

ix = colData(x)$sample_id == samples[3]
col[ix] <- col[ix] + 570 

ix = colData(x)$sample_id == samples[4]
row[ix] <- row[ix] + 79 
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
col[ix] <- col[ix] + 120

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
col[ix] <- col[ix] + 400

ix = colData(x)$sample_id == samples[2]
row[ix] <- row[ix] + 400
ix = colData(x)$sample_id == samples[2]
col[ix] <- col[ix] + 450

ix = colData(x)$sample_id == samples[3]
row[ix] <- row[ix] + 300
ix = colData(x)$sample_id == samples[3]
col[ix] <- col[ix] + 270


colData(x)$row <- row
colData(x)$col <- col

spe = cbind(spe,x)

# pal <- unname(palette.colors(12, palette = "Okabe-Ito"))
# check offsets give non-overlapping samples
df <- cbind.data.frame(colData(spe), spatialCoords(spe))
ggplot(df, aes(x = row, y = col, color = sample_id)) + 
  geom_point(size = 1) + 
  #scale_color_manual(values = cols) + 
  coord_fixed() + 
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  theme_bw()

}