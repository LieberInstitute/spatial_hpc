####################################
# spatial_HPC project
# CP markers across PRECAST clusters
# Anthony Ramnauth, Dec 19 2022
####################################

setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(here)
  library(scater)
  library(scran)
  library(sessioninfo)
})

# Load SPE
load(file = here::here("processed-data", "06_clustering", "PRECAST", "spe_modify_PRECAST_k15.Rdata"))

# Creat logcounts in assays slot
set.seed(20220201)
spe$scran_quick_cluster <- quickCluster(
  spe,
  block = spe$sample_id
)

spe <-
  computeSumFactors(spe,
                    clusters = spe$scran_quick_cluster
  )

table(spe$scran_quick_cluster)

summary(sizeFactors(spe))

spe <- logNormCounts(spe)

# set gene names as row names for easier plotting
rownames(spe) <- rowData(spe)$gene_name

pdf(file = here::here("plots", "06_clustering", "PRECAST", "CP_markers_PRECAST.pdf"))

plotExpression(spe, features=c("MSX1", "CAPS", "TTR", "KCNJ13", "TPM2", "MGP",
    "TAGLN", "SLC13A4", "ISYNA1", "FHIT", "MYL9"),
    x="PRECAST_k15", colour_by="PRECAST_k15")

dev.off()
