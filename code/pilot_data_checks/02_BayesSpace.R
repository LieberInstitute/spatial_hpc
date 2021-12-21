#tutorial from: https://edward130603.github.io/BayesSpace/articles/joint_clustering.html

setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
  library("here")
  library("spatialLIBD")
  library("ggplot2")
  library("patchwork")
  library("scater")
  library("harmony")
  library("BayesSpace")
  library("scran")
})

# load SPE
load(file=here::here("processed-data","pilot_data_checks","spe.Rdata"))
dim(spe)
# [1] 27633 28871

# Pre process
spePP = spatialPreprocess(spe, n.PCs = 50)
spePP = runUMAP(spePP, dimred = "PCA")

# Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
# Also defined by ‘spam’
# Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
# Also defined by ‘spam’
# Warning message:
#   In .check_reddim_names(x, value, withDimnames) :
#   non-NULL 'rownames(value)' should be the same as 'colnames(x)' for
# 'reducedDim<-'. This will be an error in the next release of
# Bioconductor.

colnames(reducedDim(spePP, "UMAP")) = c("UMAP1", "UMAP2")

# UMAP plots
dir.create(here::here("processed-data", "pilot_data_checks", "plots"), showWarnings = FALSE)
pdf(file=here::here("processed-data", "pilot_data_checks", "plots", "hpc_UMAP_spe.pdf"))
ggplot(data.frame(reducedDim(spePP, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spePP$sample_id))) +
  geom_point() +
  labs(color = "Sample/Capture Area") +
  theme_bw()

ggplot(data.frame(reducedDim(spePP, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spePP$subject))) +
  geom_point() +
  labs(color = "Subject/Brain") +
  theme_bw()

ggplot(data.frame(reducedDim(spePP, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(spePP$slice))) +
  geom_point() +
  labs(color = "Tissue Slices") +
  theme_bw()

dev.off()

# batch correction - HARMONY
# install.packages("devtools")
# devtools::install_github("immunogenomics/harmony")

speH = RunHarmony(spePP, "sample_id", verbose = F)
speH = runUMAP(speH, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(speH, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")

pdf(file=here::here("processed-data", "pilot_data_checks", "plots", "hpc_UMAP_harmony.pdf"))
ggplot(data.frame(reducedDim(speH, "UMAP.HARMONY")), 
       aes(x = UMAP1, y = UMAP2, color = factor(speH$sample_id))) +
  geom_point() +
  labs(color = "Sample/Capture Area") +
  theme_bw()

ggplot(data.frame(reducedDim(speH, "UMAP.HARMONY")), 
       aes(x = UMAP1, y = UMAP2, color = factor(speH$subject))) +
  geom_point() +
  labs(color = "Subject/Brain") +
  theme_bw()

ggplot(data.frame(reducedDim(speH, "UMAP.HARMONY")), 
       aes(x = UMAP1, y = UMAP2, color = factor(speH$slice))) +
  geom_point() +
  labs(color = "Tissue Slices") +
  theme_bw()

dev.off()

save(speH, file = here::here("processed-data","pilot_data_checks","spe_harmony.Rdata"))

# OFFSET
# The spatial locations for each sample need to be offset so that spots of different samples are not neighbors. 
# summary(spatialData(speH)$array_row)
# summary(spatialData(speH)$array_col)
auto_offset_row <- as.numeric(factor(unique(speH$sample_id))) * 100
names(auto_offset_row) <-unique(speH$sample_id)
speH$row <- spatialData(speH)$array_row + auto_offset_row[speH$sample_id]
speH$col <- spatialData(speH)$array_col
# summary(colData(speH)$row)
# summary(colData(speH)$col)

pdf(file=here::here("processed-data", "pilot_data_checks", "plots", "hpc_BayesSpace_OffsetCheck.pdf"))
clusterPlot(speH, "subject", color = NA) + #make sure no overlap between samples
  labs(fill = "Subject", title = "Offset check")
dev.off()

Sys.time() 
speB = spatialCluster(speH, use.dimred = "HARMONY", q = 7, nrep = 10000) #use HARMONY
Sys.time()

# [1] "2021-12-21 14:00:56 EST"

pdf(file=here::here("processed-data", "pilot_data_checks", "plots", "hpc_BayesSpace_clusterPlot.pdf"))
clusterPlot(spe, color = NA) + #plot clusters
  labs(title = "BayesSpace joint clustering")
dev.off()

save(speB, file=here::here("processed-data","pilot_data_checks", "spe_bayesSpace.Rdata"))