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
load(file = here::here("processed-data", "06_Clustering", "BayesSpace", "1st_run", "spe_modify.Rdata"), verbose = TRUE)

## load clusters
load(file = here::here("processed-data", "06_Clustering", "mbkmeans", "1st_run", "mbkmeans.Rdata"), verbose = TRUE)
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
source(file = here::here("code", "06_Clustering", "mbkmeans", "1st_run", "calculate_adj.R"))
scale.fac = 0.08139677
adj = calculate.adj.matrix(x_pixel,y_pixel,x_pixel,y_pixel,image=img,histology=TRUE)
refined = refine(sample_id = "V11U08-081_C1", pred = speb$kmeans, dis = adj, shape = "hexagon")
speb$kmeans_refine = refined
vis_clus(spe = speb, sampleid = "V11U08-081_C1", clustervar = "kmeans_refine", colors = cols, point_size = 1)

source(file = here::here("code", "06_Clustering", "mbkmeans", "1st_run", "refine.R"))
samples = unique(spe$sample_id)

for (i in 1:length(samples)){
  
  speb = spe[, which(spe$sample_id == samples[i])]
  x = spatialCoords(speb)[,1]
  y = spatialCoords(speb)[,2]
  #x = colData(speb)$array_row
  #y = colData(speb)$array_col
  
  mat = as.data.frame(cbind(x,y))
  adj = as.matrix(dist(mat, method= "euclidean", diag = TRUE, upper = TRUE))
  
  refined = refine(sample_id = colnames(speb), pred = speb$kmeans, dis = adj, shape = "hexagon")
  speb$kmeans_refine = refined
  if (i == 1) {
    sper <- speb
  } else {
    sper <- cbind(sper, speb)
  }
}

sper$annotations = "none"
sper$annotations[which(sper$kmeans_refine == "1")]="SLM"
sper$annotations[which(sper$kmeans_refine == "2")]="SLM"
sper$annotations[which(sper$kmeans_refine == "3")]="CA1"
sper$annotations[which(sper$kmeans_refine == "4")]="crap"
sper$annotations[which(sper$kmeans_refine == "5")]="SR"
sper$annotations[which(sper$kmeans_refine == "6")]="choroid"
sper$annotations[which(sper$kmeans_refine == "7")]="CA2_3"
sper$annotations[which(sper$kmeans_refine == "8")]="WM"
sper$annotations[which(sper$kmeans_refine == "9")]="choroid"
sper$annotations[which(sper$kmeans_refine == "10")]="CA1"
sper$annotations[which(sper$kmeans_refine == "11")]="ML"
sper$annotations[which(sper$kmeans_refine == "12")]="crap"
sper$annotations[which(sper$kmeans_refine == "13")]="SO"
sper$annotations[which(sper$kmeans_refine == "14")]="SLM"
sper$annotations[which(sper$kmeans_refine == "15")]="GCL"
sper$annotations[which(sper$kmeans_refine == "16")]="SUB/SL"
sper$annotations[which(sper$kmeans_refine == "17")]="CA4"

sper$kmeans_refine = paste0("mkr",sper$kmeans_refine)

spe = sper

brains = c("Br6423","Br6432","Br2743","Br8325","Br3942","Br6471","Br8667","Br8492","Br6522")
cols <- Polychrome::palette36.colors(17)
names(cols) <- sort(unique(spe$kmeans_refine))

pdf(here("plots", "06_Clustering", "mbkmeans", "1st_run", paste0("mbkmeans17_smoothed.pdf")), width = 21, height = 20)
# QC plot of tissue spots discarded

ii <- 1
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples
p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 2
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p3, p2, p4, nrow = 2)

##
ii <- 3
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

p1 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p3, p2, p4, nrow = 2)

##
ii <- 4
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[5], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p5 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
grid.arrange(p1, p2, p3, p4, p5, nrow = 2)

##
ii <- 5
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 6
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, nrow = 2)

##
ii <- 7
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 8
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 9
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans_refine", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

