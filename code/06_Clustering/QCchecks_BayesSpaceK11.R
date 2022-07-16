setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
    library(here)
    library(sessioninfo)
    library(SpatialExperiment)
    library(spatialLIBD)
    library(BayesSpace)
    library(ggplot2)
    library(Polychrome)
    library(tidySingleCellExperiment)
})
suppressPackageStartupMessages(library("ggspavis"))
suppressPackageStartupMessages(library("gridExtra"))

load(file = here::here("processed-data", "06_Clustering", "spe_modify.Rdata"))
spe <- speB
brains <- unique(spe$brnum)
# brains = c("Br6423","Br6432","Br2743","Br8325","Br3942","Br6471","Br8667","Br8492","Br6522")
pdf(here("plots", "06_Clustering", "QCchecks_BayesSpaceK11", "brains.pdf"), width = 8, height = 10)

##
ii <- 1
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = FALSE)
p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = FALSE)
p3 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = FALSE)
p4 <- plotVisium(speb[, which(speb$sample_id == samples[4])], spots = FALSE)
grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 2
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = FALSE)
p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = FALSE)

grid.arrange(p1, p2, nrow = 2)

##
ii <- 3
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = FALSE)
p2 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = FALSE)

grid.arrange(p1, p2, nrow = 1)

p1 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = FALSE)
p2 <- plotVisium(speb[, which(speb$sample_id == samples[4])], spots = FALSE)

grid.arrange(p1, p2, nrow = 2)

##
ii <- 4
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = FALSE)
p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = FALSE)
p3 <- plotVisium(speb[, which(speb$sample_id == samples[5])], spots = FALSE)
p4 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = FALSE)
p5 <- plotVisium(speb[, which(speb$sample_id == samples[4])], spots = FALSE)
grid.arrange(p1, p2, p3, p4, p5, nrow = 2)

##
ii <- 5
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = FALSE)
p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = FALSE)
p3 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = FALSE)
p4 <- plotVisium(speb[, which(speb$sample_id == samples[4])], spots = FALSE)

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 6
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = FALSE)
p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = FALSE)
p3 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = FALSE)

grid.arrange(p1, p2, p3, nrow = 2)

##
ii <- 7
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = FALSE)
p2 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = FALSE)
p3 <- plotVisium(speb[, which(speb$sample_id == samples[4])], spots = FALSE)
p4 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = FALSE)

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 8
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = FALSE)
p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = FALSE)

grid.arrange(p1, p2, nrow = 1)

##
ii <- 9
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- plotVisium(speb[, which(speb$sample_id == samples[1])], spots = FALSE)
p2 <- plotVisium(speb[, which(speb$sample_id == samples[2])], spots = FALSE)
p3 <- plotVisium(speb[, which(speb$sample_id == samples[3])], spots = FALSE)
p4 <- plotVisium(speb[, which(speb$sample_id == samples[4])], spots = FALSE)

grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()


##

k <- 11
cols <- Polychrome::palette36.colors(k)
names(cols) <- sort(unique(spe$BayesSpace_harmony_k11_nrep10000))

pdf(here("plots", "06_Clustering", "QCchecks_BayesSpaceK11", "BayesSpace_clusters_K11.pdf"), width = 21, height = 20)
# QC plot of tissue spots discarded

ii <- 1
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples
p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 2
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p3, p2, p4, nrow = 2)

##
ii <- 3
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

p1 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p3, p2, p4, nrow = 2)

##
ii <- 4
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[5], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p5 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
grid.arrange(p1, p2, p3, p4, p5, nrow = 2)

##
ii <- 5
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 6
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, nrow = 2)

##
ii <- 7
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 8
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 9
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "BayesSpace_harmony_k11_nrep10000", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

dev.off()

#########
pdf(here("plots", "06_Clustering", "QCchecks_BayesSpaceK11", "NECAB1.pdf"), width = 21, height = 20)
# QC plot of tissue spots discarded
gene <- # "NECAB1; ENSG00000123119" , "MPPED1; ENSG00000186732", "SLC17A6; ENSG00000091664" , "PROX1; ENSG00000117707"

    ii <- 1
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples
p1 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_gene(spe = speb, sampleid = samples[4], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 2
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p3, p2, p4, nrow = 2)

##
ii <- 3
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

p1 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_gene(spe = speb, sampleid = samples[4], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_gene(spe = speb, sampleid = samples[4], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p3, p2, p4, nrow = 2)

##
ii <- 4
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_gene(spe = speb, sampleid = samples[5], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p5 <- vis_gene(spe = speb, sampleid = samples[4], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
grid.arrange(p1, p2, p3, p4, p5, nrow = 2)

##
ii <- 5
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_gene(spe = speb, sampleid = samples[4], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 6
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, nrow = 2)

##
ii <- 7
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_gene(spe = speb, sampleid = samples[4], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 8
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 9
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_gene(spe = speb, sampleid = samples[4], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

dev.off()

## Cluster assignments
# 11 = crap tissue
# 9 = GCL, dendrites (molecular layer)
# 1 = white matter
# 3 = sub granular layer (inter neurons)
# 4 = CA3
# 2,5,8 different layers of CA3
# 7 = CA1/CA2
