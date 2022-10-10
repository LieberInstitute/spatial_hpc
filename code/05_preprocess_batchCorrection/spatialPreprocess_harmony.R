# tutorial from: https://edward130603.github.io/BayesSpace/articles/joint_clustering.html

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
    library("schex")
})

# load SPE
load(file = here::here("processed-data", "04_QC", "spe_QC.Rdata"))
dim(spe)
# [1] 30359 135640

spe$discard_auto_br <- spe$low_sum_br | spe$low_detected_br
spe <- spe[, colData(spe)$discard_auto_id == FALSE]
dim(spe)

# Pre process
spe <- spatialPreprocess(spe, n.PCs = 50)

message("Running runUMAP()")
Sys.time()
set.seed(20220208)
spe <- runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) <- c("UMAP1", "UMAP2")
Sys.time()
# Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
# Also defined by ‘spam’
# Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
# Also defined by ‘spam’
# Warning message:
#   In .check_reddim_names(x, value, withDimnames) :
#   non-NULL 'rownames(value)' should be the same as 'colnames(x)' for
# 'reducedDim<-'. This will be an error in the next release of
# Bioconductor.


hex <- make_hexbin(spe, nbins = 100, dimension_reduction = "UMAP", use_dims = c(1, 2))

# UMAP plots
pdf(file = here::here("plots", "05_Batch_correction", "spatialPreprocess_UMAP.pdf"))
ggplot(
    data.frame(reducedDim(spe, "UMAP")),
    aes(x = UMAP1, y = UMAP2, color = factor(spePP$brnum))) +
    geom_point() +
    labs(color = "Subject/Brain") +
    theme_bw()

label_df <- make_hexbin_label(hex, col = "brnum")
plot_hexbin_meta(hex, col = "brnum", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("Brains") + theme(legend.position = "right")

ggplot(
    data.frame(reducedDim(spe, "UMAP")),
    aes(x = UMAP1, y = UMAP2, color = factor(spePP$sample_id))) +
    geom_point() +
    labs(color = "Sample/Capture Area") +
    theme_bw()

label_df <- make_hexbin_label(hex, col = "sample_id")
plot_hexbin_meta(hex, col = "sample_id", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("Capture area") + theme(legend.position = "right")

dev.off()

# batch correction - HARMONY
# install.packages("devtools")
# devtools::install_github("immunogenomics/harmony")

spe <- RunHarmony(spe, "sample_id", verbose = F)
spe <- runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(spe, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")

hex <- make_hexbin(spe, nbins = 100, dimension_reduction = "UMAP.HARMONY", use_dims = c(1, 2))

pdf(file = here::here("plots", "05_Batch_correction", "spatialPreprocess_UMAP_harmony.pdf"))
ggplot(
    data.frame(reducedDim(spe, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(speH$brnum))) +
    geom_point() +
    labs(color = "Subject/Brain") +
    theme_bw()

label_df <- make_hexbin_label(hex, col = "brnum")
plot_hexbin_meta(hex, col = "brnum", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("Brains") + theme(legend.position = "right")

ggplot(
    data.frame(reducedDim(spe, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(speH$sample_id))) +
    geom_point() +
    labs(color = "Sample/Capture Area") +
    theme_bw()

label_df <- make_hexbin_label(hex, col = "sample_id")
plot_hexbin_meta(hex, col = "sample_id", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("Capture area") + theme(legend.position = "right")

dev.off()

save(spe, file = here::here("processed-data", "05_Batch_correction", "spatialPreprocess_spe.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()