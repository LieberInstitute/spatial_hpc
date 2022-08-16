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
spe$kmeans <- km_res[[17]]$Clusters

colData(spe)$imagerow = spatialCoords(spe)[,2]
colData(spe)$imagecol = spatialCoords(spe)[,1]

spec = spe[, which(spe$brnum == "Br3942")]
spec = spatialEnhance(spec, q = 17, platform = "Visium", use.dimred = "HARMONY", nrep = 100, init = "kmeans", burn.in = 10)
