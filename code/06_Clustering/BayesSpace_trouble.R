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
suppressPackageStartupMessages(library("gridExtra"))

# Load SPE
load(file = here::here("processed-data", "06_Clustering", "spe_modify.Rdata"))
unique(spe$brnum)

load(file = here("processed-data", "06_Clustering", "mbkmeans.Rdata"))
spe$mbkmeans <- km_res[[13]]$Clusters

speo = spe
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

# df = cbind.data.frame(colData(x))
# ggplot(df, aes(x = array_row, y = array_col)) + geom_point(size = 1) 

speb = speb[,colData(speb)$sample_id != "V11L05-335_D1"]
speb = cbind(speb,x)

# The spatial locations for each sample need to be offset so that spots of different samples are not neighbors.
# summary(colData(spe)$array_row)
# summary(colData(spe)$array_col)

row <- colData(speb)$array_row
col <- colData(speb)$array_col

ix = colData(speb)$sample_id == "V11A20-297_B1"
row[ix] <- row[ix] + 75 

ix = colData(speb)$sample_id == "V11A20-297_C1"
row[ix] <- row[ix] - 10 

ix = colData(speb)$sample_id == "V11A20-297_C1"
col[ix] <- col[ix] + 130 

ix = colData(speb)$sample_id == "V11L05-335_D1"
row[ix] <- row[ix] + 70 

ix = colData(speb)$sample_id == "V11L05-335_D1"
col[ix] <- col[ix] + 130 

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
k <- 17

# Run BayesSpace
message("Running spatialCluster()")
Sys.time()
set.seed(12345)
spe <- spatialCluster(speb, use.dimred = "HARMONY", q = k, platform = "Visium", nrep = 50000)
Sys.time()

save(spe, file = here("processed-data", "06_Clustering", "BayesSpace_rerun_trouble_k17.Rdata"))

samples <- unique(spe$sample_id)
samples

cols <- Polychrome::palette36.colors(k)


pdf(here("plots", "06_Clustering", "BayesSpace", "BayesSpace_rerun_trouble_k17.pdf"), width = 21, height = 20)
clustV = "cluster.init"
names(cols) <- sort(unique(spe$cluster.init))
p1 <- vis_clus(spe = spe, sampleid = samples[1], clustervar = clustV, colors = cols, point_size = 2)
p2 <- vis_clus(spe = spe, sampleid = samples[4], clustervar = clustV, colors = cols, point_size = 2)
p3 <- vis_clus(spe = spe, sampleid = samples[2], clustervar = clustV, colors = cols, point_size = 2)
p4 <- vis_clus(spe = spe, sampleid = samples[3], clustervar = clustV, colors = cols, point_size = 2)
grid.arrange(p1, p2, p3, p4, nrow = 2)

clustV = "spatial.cluster"
names(cols) <- sort(unique(spe$spatial.cluster))
p1 <- vis_clus(spe = spe, sampleid = samples[1], clustervar = clustV, colors = cols, point_size = 2)
p2 <- vis_clus(spe = spe, sampleid = samples[4], clustervar = clustV, colors = cols, point_size = 2)
p3 <- vis_clus(spe = spe, sampleid = samples[2], clustervar = clustV, colors = cols, point_size = 2)
p4 <- vis_clus(spe = spe, sampleid = samples[3], clustervar = clustV, colors = cols, point_size = 2)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
