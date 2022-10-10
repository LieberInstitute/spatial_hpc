## Load remaining required packages
setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

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


## Create output directories
## Load the data
load(file = here::here("processed-data", "04_QC", "spe_QC.Rdata"))
dim(spe)

spe$discard_auto_br <- spe$low_sum_br | spe$low_detected_br
spe <- spe[, colData(spe)$discard_auto_id == FALSE]
dim(spe)

message("Running quickCluster()")
set.seed(20220201)
Sys.time()
spe$scran_quick_cluster <- quickCluster(
  spe,
  BPPARAM = MulticoreParam(4),
  block = spe$sample_id,
  block.BPPARAM = MulticoreParam(4)
)
Sys.time()

message("Running computeSumFactors()")
Sys.time()
## Might be needed:
# options(error = recover)
spe <-
  computeSumFactors(spe,
                    clusters = spe$scran_quick_cluster,
                    BPPARAM = MulticoreParam(4)
  )
Sys.time()

table(spe$scran_quick_cluster)

message("Running checking sizeFactors()")
summary(sizeFactors(spe))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.006964  0.368746  0.662599  1.000000  1.223521 28.347141 

message("Running logNormCounts()")
spe <- logNormCounts(spe)


message("Running modelGeneVar()")
## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
dec <- modelGeneVar(spe,
                    block = spe$sample_id,
                    BPPARAM = MulticoreParam(4)
)

pdf(file = here::here("plots", "05_Batch_correction", "scran_modelGeneVar_captureArea.pdf"), useDingbats = FALSE)
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
#[1] 1804

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
#[1] 14571

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
#[1] 13820

save(top.hvgs,
     top.hvgs.fdr5,
     top.hvgs.fdr1,
     file = here::here("processed-data", "05_Batch_correction", "top.hvgs_captureArea.Rdata")
)


message("Running runPCA()")
set.seed(20220201)
Sys.time()
spe <- runPCA(spe, subset_row = top.hvgs, ncomponents = 50)
Sys.time()

# make elbow plot to determine PCs to use
percent.var <- attr(reducedDim(spe, "PCA"), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow

pdf(file = here::here("plots", "05_Batch_correction", "pca_elbow_captureArea.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()

# message("Running runTSNE() perplexity 5") ignore takes long time
# Sys.time()
# set.seed(20220201)
# spe <-
#   runTSNE(spe,
#           dimred = "PCA",
#           name = "TSNE_perplexity05",
#           perplexity = 5
#   )
# Sys.time()


message("Running runUMAP()")
Sys.time()
set.seed(20220208)
spe <- runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) <- c("UMAP1", "UMAP2")
Sys.time()

## Perform harmony batch correction
message("Running RunHarmony()")
Sys.time()
set.seed(20220208)
spe <- RunHarmony(spe, "sample_id", verbose = FALSE)
Sys.time()

message("Running runUMAP() on HARMONY dimensions")
Sys.time()
set.seed(20220208)
spe <- runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")
Sys.time()
colnames(reducedDim(spe, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")


# message("Running runTSNE() perplexity 5 on HARMONY dimensions")
# Sys.time()
# set.seed(20220208)
# spe <-
#   runTSNE(spe,
#           dimred = "HARMONY",
#           name = "TSNE_perplexity05.HARMONY",
#           perplexity = 5
#   )
# Sys.time()


## Explore UMAP results
pdf(file = here::here("plots", "05_Batch_correction", "OSCAPreprocess_captureArea_UMAP.pdf"))
ggplot(
  data.frame(reducedDim(spe, "UMAP")),
  aes(x = UMAP1, y = UMAP2, color = factor(spe$brnum))
) +
  geom_point() +
  labs(color = "brnum") +
  theme_bw()

hex <- make_hexbin(spe, nbins = 100, dimension_reduction = "UMAP", use_dims = c(1, 2))
label_df <- make_hexbin_label(hex, col = "brnum")
plot_hexbin_meta(hex, col = "brnum", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("Brain") + theme(legend.position = "right")

ggplot(
  data.frame(reducedDim(spe, "UMAP")),
  aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))
) +
  geom_point() +
  labs(color = "sample_id") +
  theme_bw()

label_df <- make_hexbin_label(hex, col = "sample_id")
plot_hexbin_meta(hex, col = "sample_id", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("Capture area") + theme(legend.position = "right")

dev.off()

## Explore UMAP on HARMONY reduced dimensions
pdf(file = here::here("plots", "05_Batch_correction", "OSCAPreprocess_captureArea_UMAP_harmony.pdf"), width = 9)
ggplot(
  data.frame(reducedDim(spe, "UMAP.HARMONY")),
  aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))
) +
  geom_point() +
  labs(color = "brnum") +
  theme_bw()

hex <- make_hexbin(spe, nbins = 100, dimension_reduction = "UMAP.HARMONY", use_dims = c(1, 2))
label_df <- make_hexbin_label(hex, col = "brnum")
plot_hexbin_meta(hex, col = "brnum", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("HARMONY Brains") + theme(legend.position = "right")

ggplot(
  data.frame(reducedDim(spe, "UMAP.HARMONY")),
  aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))
) +
  geom_point() +
  labs(color = "sample_id") +
  theme_bw()

label_df <- make_hexbin_label(hex, col = "sample_id")
plot_hexbin_meta(hex, col = "sample_id", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("HARMONY Capture area") + theme(legend.position = "right")
dev.off()

save(spe, file = here::here("processed-data", "05_Batch_correction", "OSCAPreprocess_captureArea_spe.Rdata"))


## Object size in GB
## (do this near the end in case lobstr crashes, it's happened to me once)
lobstr::obj_size(spe)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()