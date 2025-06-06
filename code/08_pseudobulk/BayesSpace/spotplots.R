setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
  library(here)
  library(sessioninfo)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(BayesSpace)
  library(ggplot2)
})
suppressPackageStartupMessages(library("ggspavis"))
suppressPackageStartupMessages(library("gridExtra"))

load(file = here::here("processed-data", "06_Clustering", "spe_modify.Rdata"))

brains = c("Br6423","Br6432","Br2743","Br8325","Br3942","Br6471","Br8667","Br8492","Br6522")
SVGs <- c(
  "TNNT2",
  "GRM2",
  "POSTN",
  "PDYN",
  "GRP",
  "MCOLN3",
  "SLC30A3",
  "KCNIP2",
  "P2RX2",
  "AC106795.2"
)
SVGs <- c(
  "CDKN1C",
  "RNASE1",
  "BCAS1",
  "ENPP2")
# Locate the marker genes
SVG_search <- rowData(spe)$gene_search[match(SVGs, rowData(spe)$gene_name)]

for (i in SVG_search) {
  gene <- i
  gene_name <- strsplit(i, ";")[[1]][1]
  pdf(here("plots", "08_pseudobulk", "BayesSpace", "DG",paste0(gene_name, ".pdf")), width = 21, height = 20)
  
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
}
