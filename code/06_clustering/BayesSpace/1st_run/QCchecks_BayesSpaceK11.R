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

load(file = here::here("processed-data", "06_Clustering", "BayesSpace", "1st_run", "spe_modify.Rdata"))

brains <- unique(spe$brnum)
# brains = c("Br6423","Br6432","Br2743","Br8325","Br3942","Br6471","Br8667","Br8492","Br6522")
pdf(here("plots", "06_Clustering", "BayesSpace", "1st_run", "QCchecks_BayesSpaceK11", "brains.pdf"), width = 8, height = 10)

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

pdf(here("plots", "06_Clustering", "BayesSpace", "1st_run", "QCchecks_BayesSpaceK11", "BayesSpace_clusters_K11.pdf"), width = 21, height = 20)
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
pdf(here("plots", "06_Clustering", "BayesSpace", "1st_run", "QCchecks_BayesSpaceK11", "NECAB1.pdf"), width = 21, height = 20)
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

table(spe$BayesSpace_harmony_k11_nrep10000, spe$brnum)

#     Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325 Br8492 Br8667
# 1    1683   2064   1743   1094   1636   3034   2499   1127   2005
# 2     781   1957   1270    724    750   1166   2431    335   1244
# 3      45   1179     61     65   1560    245    879   1401   3202
# 4    4601   1741   5379   1582   1682   3071   2959    109   1006
# 5    1693   2238   1673   1269   2754   3602   2536   1206   2870
# 6    1277   2388   1078    888   1279   2471   1663    688   1316
# 7    1499   2604    806   1412    726    945   1863    841   2820
# 8    1983   3473   1137   1124   1530   1920   3492    703   1409
# 9      68    817     66    226   1304    448    647   1203    986
# 10    282    590    479    150    739    549    716    568    468
# 11    175    389    311     85    271   1189    730    191    507
