
setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

library("SpatialExperiment")
library("here")
library("sessioninfo")
library("jaffelab")
library("ggplot2")
library("spatialLIBD")
library("Polychrome")
suppressPackageStartupMessages(library("gridExtra"))

load(file = here::here("processed-data", "06_clustering", "BayesSpace", "1st_run",
    "spe_modify.Rdata"))

brains = c("Br6423","Br6432","Br2743","Br8325","Br3942","Br6471","Br8667","Br8492","Br6522")
cols <- Polychrome::palette36.colors(16)
names(cols) <- sort(unique(spe$ManualAnnotation))

pdf(here("plots", "manual_annotation", "manual_annotation_test.pdf"), width = 21, height = 20)
# QC plot of tissue spots discarded

ii <- 1
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples
p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 2
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p3, p2, p4, nrow = 2)

##
ii <- 3
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

p1 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p3, p2, p4, nrow = 2)

##
ii <- 4
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[5], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p5 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
grid.arrange(p1, p2, p3, p4, p5, nrow = 2)

##
ii <- 5
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 6
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, nrow = 2)

##
ii <- 7
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 8
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 9
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "ManualAnnotation", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

dev.off()


