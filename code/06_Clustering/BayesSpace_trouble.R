###############################
# Set up SGE array job to run k=5 to k = 20
# Found in BayesSpaces.sh shell script line -t 2-15
setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
  library(here)
  library(sessioninfo)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(BayesSpace)
  library(ggplot2)
  library(Polychrome)
})

# Load SPE
load(file = here::here("processed-data", "06_Clustering", "spe_modify.Rdata"))
unique(spe$brnum)

speb <- spe[, which(spe$brnum == "Br8325")]
unique(speb$sample_id)

speb = speb[,which(speb$sample_id != "V11A20-297_D1")]

# rotate array coordinates

radians <- degrees * pi / 180
refl_vec <- c(1, 1)

rotation_mat <- matrix(c(cos(radians), sin(radians), -1 * sin(radians), cos(radians)),nrow = 2)

dim_max <- dim(imgRaster(x)) / scaleFactors(x)[1]
new_coords <- refl_vec * t(spatialCoords(x)) 
new_coords <- rotation_mat %*% new_coords
new_coords <- new_coords + rev(dim_max)
new_coords <- t(new_coords + c(0, 0)) # transpose and translate
colnames(new_coords) <- colnames(spatialCoords(x))

#   Ensure points are at integer values
new_coords <- round(new_coords)


# Choose k
k <- 15

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(speb)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

# The spatial locations for each sample need to be offset so that spots of different samples are not neighbors.
# summary(colData(spe)$array_row)
# summary(colData(spe)$array_col)
sample2 = "V11A20-297_C1"
sample2 = "V11L05-335_D1"
sample1 = "V11A20-297_A1"
sample4 = "V11A20-297_B1"

row <- colData(speb)$array_row
col <- colData(speb)$array_col

ix = colData(speb)$sample_id == "V11A20-297_B1"
row[ix] <- row[ix] + 100 

ix = colData(speb)$sample_id == "V11A20-297_C1"
col[ix] <- col[ix] + 150 

ix = colData(speb)$sample_id == "V11L05-335_D1"
row[ix] <- row[ix] + 100 

ix = colData(speb)$sample_id == "V11L05-335_D1"
col[ix] <- col[ix] + 150 

pal <- unname(palette.colors(4, palette = "Okabe-Ito"))

colData(speb)$row <- row
colData(speb)$col <- col

# check offsets give non-overlapping samples
df <- cbind.data.frame(colData(speb), spatialCoords(speb))
ggplot(df, aes(x = row, y = col, color = sample_id)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = pal) + 
  coord_fixed() + 
  ggtitle("BayesSpace: multiple samples offsets check") + 
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  theme_bw()

# Run BayesSpace
message("Running spatialCluster()")
Sys.time()
set.seed(12345)
spe <- spatialCluster(spe, use.dimred = "HARMONY", q = k, platform = "Visium", save.chain = TRUE, nrep = 10000)
Sys.time()

save(spe, file = here("processed-data", "06_Clustering", "BayesSpace_rerun_k15.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
