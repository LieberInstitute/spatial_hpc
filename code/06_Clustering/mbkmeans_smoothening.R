setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

library("SpatialExperiment")
library("spatialLIBD")
library("mbkmeans")
library("here")
library("sessioninfo")
library("jaffelab")
library("ggplot2")
library("Polychrome")
library("BayesSpace")
suppressPackageStartupMessages(library("gridExtra"))

## load data
load(file = here::here("processed-data", "06_Clustering","spe_modify.Rdata"), verbose = TRUE)

## load clusters
load(file = here::here("processed-data", "06_Clustering","mbkmeans.Rdata"), verbose = TRUE)
spe$kmeans <- km_res[[13]]$Clusters

colData(spe)$imagerow = spatialCoords(spe)[,2]
colData(spe)$imagecol = spatialCoords(spe)[,1]

spec = spe[, which(spe$brnum == "Br3942")]
spec = spatialEnhance(spec, q = 17, platform = "Visium", use.dimred = "HARMONY", nrep = 100, init = "kmeans", burn.in = 10)

speb = spe[, which(spe$sample_id == "V11U08-081_C1")]
x = spatialCoords(speb)[,1]
y = spatialCoords(speb)[,2]
mat = as.data.frame(cbind(x,y))
adj = as.matrix(dist(mat, method= "euclidean", diag = TRUE, upper = TRUE))
refined = refine(sample_id = "V11U08-081_C1", pred = speb$kmeans, dis, shape = "hexagon")
speb$kmeans_refine = refined

vis_clus(spe = speb, sampleid = "V11U08-081_C1", clustervar = "kmeans_refine", colors = cols, point_size = 2)
