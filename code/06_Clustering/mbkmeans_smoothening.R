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
x = colData(speb)$array_row
y = colData(speb)$array_col

mat = as.data.frame(cbind(x,y))
adj = as.matrix(dist(mat, method= "euclidean", diag = TRUE, upper = TRUE))

refined = refine(sample_id = colnames(speb), pred = speb$kmeans, dis = adj, shape = "hexagon")
speb$kmeans_refine = refined
vis_clus(spe = speb, sampleid = "V11U08-081_C1", clustervar = "kmeans_refine", point_size = 2)

pairwise.distance <- function(X) {
  n <- nrow(X)
  adj <- matrix(rep(0, n*n),nrow = n, ncol=n)
  for (i in 1:n) {
    for (j in 1:n){
      adj[i,j] <- sqrt(sum((X[i,]-X[j,])^2))
    }
  }
  return(adj)
}

adj = pairwise.distance(mat)
refined = refine(sample_id = "V11U08-081_C1", pred = speb$kmeans, dis = adj, shape = "hexagon")
speb$kmeans_refine = refined
vis_clus(spe = speb, sampleid = "V11U08-081_C1", clustervar = "kmeans_refine", colors = cols, point_size = 2)

x = speb$array_col
y = speb$array_row
x_pixel = spatialCoords(speb)[,1]*0.08139677
y_pixel = spatialCoords(speb)[,2]*0.08139677
img = imgRaster(speb, sample_id = "V11U08-081_C1", image_id = "hires")
source(file = here::here("code", "06_Clustering","calculate_adj.R"))
scale.fac = 0.08139677
adj = calculate.adj.matrix(x_pixel,y_pixel,x_pixel,y_pixel,image=img,histology=TRUE)
refined = refine(sample_id = "V11U08-081_C1", pred = speb$kmeans, dis = adj, shape = "hexagon")
speb$kmeans_refine = refined
vis_clus(spe = speb, sampleid = "V11U08-081_C1", clustervar = "kmeans_refine", colors = cols, point_size = 1)

source(file = here::here("code", "06_Clustering","refine.R"))