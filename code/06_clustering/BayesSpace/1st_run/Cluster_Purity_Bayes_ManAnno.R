############################################################
# spatial_HPC project
# Cluster Purity check for BayesSpace & Manual Annotations
# Anthony Ramnauth, Feb 22 2022
############################################################

setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(spatialLIBD)
  library(here)
  library(scuttle)
  library(scater)
  library(scran)
  library(dplyr)
  library(PCAtools)
  library(schex)
  library(viridis)
  library(sessioninfo)
  library(gridExtra)
  library(bluster)
})

# Load SPE
load(file = here::here("processed-data", "06_clustering", "BayesSpace", "1st_run",
    "BayesSpace_rerun_k15.Rdata"))

# Change <NA> to character in Manual Annotations

Man_Ann <- as.vector(spe$ManualAnnotation)

Man_Ann <- replace(Man_Ann, is.na(Man_Ann), "None")

spe$ManualAnnotation <- Man_Ann
spe$ManualAnnotation <- factor(spe$ManualAnnotation)

# Check cluster purity for Manual Annotations

pure_spe_man <- neighborPurity(reducedDim(spe, "PCA"), spe$ManualAnnotation)
pure_spe_man

pure_data_man <- as.data.frame(pure_spe_man)
pure_data_man$maximum <- factor(pure_data_man$maximum)
pure_data_man$cluster <- spe$ManualAnnotation

pdf(file = here::here("plots","manual_annotation", "Cluster_Purity_ManualAnnotations.pdf"), width = 10, height = 8)

ggplot(pure_data_man, aes(x=cluster, y=purity, colour=maximum)) +
    ggbeeswarm::geom_quasirandom(method="smiley")

boxplot(split(pure_data_man$purity, pure_data_man$cluster))

dev.off()

# Check cluster purity for BayesSpace k = 15

pure_spe_bayes <- neighborPurity(reducedDim(spe, "PCA"), spe$BayesSpace_harmony_k15_nrep10000)
pure_spe_bayes

pure_data_bayes <- as.data.frame(pure_spe_bayes)
pure_data_bayes$maximum <- factor(pure_data_bayes$maximum)
pure_data_bayes$cluster <- factor(spe$BayesSpace_harmony_k15_nrep10000)

pdf(file = here::here("plots","06_clustering", "BayesSpace", "1st_run", "Cluster_Purity_BayesSpace_k15.pdf"),
    width = 10, height = 8)

ggplot(pure_data_bayes, aes(x=cluster, y=purity, colour=maximum)) +
    ggbeeswarm::geom_quasirandom(method="smiley")

boxplot(split(pure_data_bayes$purity, pure_data_bayes$cluster))

dev.off()



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
