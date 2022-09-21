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

speb <- spe[,colData(spe)$brnum == "Br8325"]
unique(speb$sample_id)

speb = speb[,colData(speb)$sample_id != "V11A20-297_D1"]

# rotate array coordinates

radians <- 180 * pi / 180
refl_vec <- c(1, 1)

rotation_mat <- matrix(c(cos(radians), sin(radians), -1 * sin(radians), cos(radians)),nrow = 2)

x = speb[,colData(speb)$sample_id == "V11L05-335_D1"]

dim_max <- c(max(x$array_row), max(x$array_col))
arrayCoords <- cbind(x$array_col, x$array_row)
rownames(arrayCoords) <- rownames(spatialCoords(x))
colnames(arrayCoords) <- c("array_col", "array_row")
new_coords <- refl_vec * t(arrayCoords) 
new_coords <- rotation_mat %*% new_coords
new_coords <- new_coords + rev(dim_max)
new_coords <- t(new_coords + c(0, 0)) # transpose and translate
colnames(new_coords) <- colnames(arrayCoords)

#   Ensure points are at integer values
new_coords <- round(new_coords)
x$array_col = new_coords[,1]
x$array_row = new_coords[,2]

ggplot(x, aes(x = x$array_row, y = x$array_col)) + geom_point(size = 1) 

speb[,colData(speb)$sample_id == "V11L05-335_D1"] = x

# The spatial locations for each sample need to be offset so that spots of different samples are not neighbors.
# summary(colData(spe)$array_row)
# summary(colData(spe)$array_col)

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

## Set the BayesSpace metadata using code from
## https://github.com/edward130603/BayesSpace/blob/master/R/spatialPreprocess.R#L43-L46
metadata(speb)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

# Choose k
k <- 15

# Run BayesSpace
message("Running spatialCluster()")
Sys.time()
set.seed(12345)
spe <- spatialCluster(speb, use.dimred = "HARMONY", q = k, platform = "Visium", nrep = 50000)
Sys.time()

save(spe, file = here("processed-data", "06_Clustering", "BayesSpace_rerun_trouble_k15.Rdata"))

samples <- unique(spe$sample_id)
samples
clustV = "cluster.init"
cols <- Polychrome::palette36.colors(k)
p1 <- vis_clus(spe = spe, sampleid = samples[1], clustervar = clustV, colors = cols, point_size = 2)
p2 <- vis_clus(spe = spe, sampleid = samples[4], clustervar = clustV, colors = cols, point_size = 2)
p3 <- vis_clus(spe = spe, sampleid = samples[2], clustervar = clustV, colors = cols, point_size = 2)
p4 <- vis_clus(spe = spe, sampleid = samples[3], clustervar = clustV, colors = cols, point_size = 2)
grid.arrange(p1, p3, p4, p5, nrow = 2)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
