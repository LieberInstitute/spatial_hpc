## Load remaining required packages
setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
  library("here")
  library("SpatialExperiment")
  library("spatialLIBD")
  library("sessioninfo")
  library("scran") ## requires uwot for UMAP
  library("uwot")
  library("scater")
  library("BiocParallel")
  library("PCAtools")
  library("ggplot2")
  library("Polychrome")
  library("harmony")
  library("schex")
})


## Load the data
load(file = here::here("processed-data", "04_QC", "spe_QC_allSamples.Rdata"))
dim(spe)
#[1]  31483 191136

spe <- spe[, colData(spe)$discard_auto_id == FALSE]
dim(spe)
#[1]  31483 188762
##Compute library factors (OSCA)
set.seed(25323)
spe <- computeLibraryFactors(spe)
summary(sizeFactors(spe))

hist(sizeFactors(spe), breaks = 1000,xlim=c(0,10))


message("Running logNormCounts()")
spe <- logNormCounts(spe)

message("Running modelGeneVar()")
## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
dec <- modelGeneVar(spe,
                    block = spe$sample_id,
                    BPPARAM = MulticoreParam(4)
)

pdf(file = here::here("plots", "05_preprocess_batchCorrection", "scran_modelGeneVar_captureArea_allSamples.pdf"), useDingbats = FALSE)
mapply(function(block, blockname) {
  plot(
    block$mean,
    block$total,
    xlab = "Mean log-expression",
    ylab = "Variance",
    main = blockname
  )
  curve(metadata(block)$trend(x),
        col = "blue",
        add = TRUE
  )
}, dec$per.block, names(dec$per.block))
dev.off()

message("Running getTopHVGs()")
top.hvgs <- getTopHVGs(dec, prop = 0.1)
length(top.hvgs)

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.00001)
length(top.hvgs.fdr1)

save(top.hvgs,
     top.hvgs.fdr5,
     top.hvgs.fdr1,
     file = here::here("processed-data", "05_preprocess_batchCorrection", "top.hvgs_captureArea_allSamples.Rdata")
)

# 
# message("Running runPCA()")
# set.seed(20220201)
# Sys.time()
# spe <- runPCA(spe, subset_row = top.hvgs, ncomponents = 50)
# Sys.time()
# 
# # make elbow plot to determine PCs to use
# percent.var <- attr(reducedDim(spe, "PCA"), "percentVar")
# chosen.elbow <- PCAtools::findElbowPoint(percent.var)
# chosen.elbow
# 
# pdf(file = here::here("plots", "05_preprocess_batchCorrection", "pca_elbow_captureArea_allSamples.pdf"), useDingbats = FALSE)
# plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
# abline(v = chosen.elbow, col = "red")
# dev.off()
# 
# # message("Running runTSNE() perplexity 5") ignore takes long time
# # Sys.time()
# # set.seed(20220201)
# # spe <-
# #   runTSNE(spe,
# #           dimred = "PCA",
# #           name = "TSNE_perplexity05",
# #           perplexity = 5
# #   )
# # Sys.time()
# 
# 
# message("Running runUMAP()")
# Sys.time()
# set.seed(20220208)
# spe <- runUMAP(spe, dimred = "PCA")
# colnames(reducedDim(spe, "UMAP")) <- c("UMAP1", "UMAP2")
# Sys.time()
# 
# ## Perform harmony batch correction
# message("Running RunHarmony()")
# Sys.time()
# set.seed(20220208)
# spe <- RunHarmony(spe, "sample_id", verbose = FALSE)
# Sys.time()
# 
# message("Running runUMAP() on HARMONY dimensions")
# Sys.time()
# set.seed(20220208)
# spe <- runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")
# Sys.time()
# colnames(reducedDim(spe, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")
# 
# 
# # message("Running runTSNE() perplexity 5 on HARMONY dimensions")
# # Sys.time()
# # set.seed(20220208)
# # spe <-
# #   runTSNE(spe,
# #           dimred = "HARMONY",
# #           name = "TSNE_perplexity05.HARMONY",
# #           perplexity = 5
# #   )
# # Sys.time()
# 
# 
# ## Explore UMAP results
# pdf(file = here::here("plots", "05_preprocess_batchCorrection", "OSCApreprocess_captureArea_UMAP_allSamples.pdf"))
# ggplot(
#   data.frame(reducedDim(spe, "UMAP")),
#   aes(x = UMAP1, y = UMAP2, color = factor(spe$brnum))
# ) +
#   geom_point() +
#   labs(color = "brnum") +
#   theme_bw()
# 
# hex <- make_hexbin(spe, nbins = 100, dimension_reduction = "UMAP", use_dims = c(1, 2))
# label_df <- make_hexbin_label(hex, col = "brnum")
# plot_hexbin_meta(hex, col = "brnum", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("Brain") + theme(legend.position = "right")
# 
# ggplot(
#   data.frame(reducedDim(spe, "UMAP")),
#   aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))
# ) +
#   geom_point() +
#   labs(color = "sample_id") +
#   theme_bw()
# 
# label_df <- make_hexbin_label(hex, col = "sample_id")
# plot_hexbin_meta(hex, col = "sample_id", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("Capture area") + theme(legend.position = "right")
# 
# dev.off()
# 
# ## Explore UMAP on HARMONY reduced dimensions
# pdf(file = here::here("plots", "05_preprocess_batchCorrection", "OSCApreprocess_captureArea_UMAP_harmony_allSamples.pdf"), width = 9)
# ggplot(
#   data.frame(reducedDim(spe, "UMAP.HARMONY")),
#   aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))
# ) +
#   geom_point() +
#   labs(color = "brnum") +
#   theme_bw()
# 
# hex <- make_hexbin(spe, nbins = 100, dimension_reduction = "UMAP.HARMONY", use_dims = c(1, 2))
# label_df <- make_hexbin_label(hex, col = "brnum")
# plot_hexbin_meta(hex, col = "brnum", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("HARMONY Brains") + theme(legend.position = "right")
# 
# ggplot(
#   data.frame(reducedDim(spe, "UMAP.HARMONY")),
#   aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))
# ) +
#   geom_point() +
#   labs(color = "sample_id") +
#   theme_bw()
# 
# label_df <- make_hexbin_label(hex, col = "sample_id")
# plot_hexbin_meta(hex, col = "sample_id", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("HARMONY Capture area") + theme(legend.position = "right")
# dev.off()
# 
# save(spe, file = here::here("processed-data", "05_preprocess_batchCorrection", "OSCApreprocess_harmony_captureArea_spe_allSamples.Rdata"))
# 
# 
# ## Object size in GB
# ## (do this near the end in case lobstr crashes, it's happened to me once)
# lobstr::obj_size(spe)
save(spe,file=here::here('processed-data','05_preprocess_batchCorrection','spe_norm.rda'
))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()