#######################################
# spatial_hpc project
# Saving PRECAST clusters to spe object
# Anthony Ramnauth, Dec 13 2022
#######################################

setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(Seurat)
    library(SpatialExperiment)
    library(PRECAST)
    library(spatialLIBD)
    library(ggplot2)
    library(gridExtra)
    library("here")
})

# Load the PRECAST objects (looks like you can only do one at a time)

load(file = here("processed-data", "06_clustering", "PRECAST", "PRECASTObj.Rdata"))

# Load the spe object

load(file = here::here("processed-data", "04_QC", "spe_QC.Rdata"))

# Leaving this bit from Maddy's script, not sure I need it.
resList <- PRECASTObj@resList
PRECASTObj <- selectModel(PRECASTObj)

# Convert to seurat object for easier format to explore

seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")
seuInt

dim(spe)

dim(seuInt@meta.data)

# Find number of barcodes that are different from spe

length(setdiff(colnames(spe), rownames(seuInt@meta.data)))

# Create temp df for the clusters from PRECAST k_15

k_15 <- data.frame(
    cluster = seuInt@meta.data$cluster,
    barcodes = rownames(seuInt@meta.data)
)

rownames(k_15) <- rownames(seuInt@meta.data)

# Create df of NAs to fill missing barcodes & in an reorder barcodes

emptyNaDF <- data.frame(matrix(NA,nrow = 712, ncol = 2))

rownames(emptyNaDF) <- setdiff(colnames(spe), rownames(seuInt@meta.data))

colnames(emptyNaDF) <- c("cluster", "barcodes")

k_15 <- rbind(k_15, emptyNaDF)

# Reorder to match colnames(spe)
k_15 <- k_15[order(match(rownames(k_15), colnames(spe))), ]

stopifnot(rownames(k_15) == colnames(spe))

# Add column of PRECAST cluster assignments to colData(spe)

colData(spe)$PRECAST_k15 <- k_15$cluster

cols <- Polychrome::palette36.colors(16)
names(cols) <- sort(unique(spe$PRECAST_k15))
brains = c("Br6423","Br6432","Br2743","Br8325","Br3942","Br6471","Br8667","Br8492","Br6522")

# Do a quick check to make sure clusters look accurately mapped to tissue

pdf(file = here("plots", "06_clustering", "PRECAST", "Test_PRECAST_k15.pdf"),
    width = 21, height = 20)

ii <- 1
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples
p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 2
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p3, p2, p4, nrow = 2)

##
ii <- 3
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

p1 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p3, p2, p4, nrow = 2)

##
ii <- 4
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[5], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p5 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
grid.arrange(p1, p2, p3, p4, p5, nrow = 2)

##
ii <- 5
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 6
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, nrow = 2)

##
ii <- 7
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 8
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 9
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "PRECAST_k15", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

dev.off()

save(spe, file = here::here("processed-data", "06_clustering", "PRECAST", "spe_modify_PRECAST_k15.Rdata"))
