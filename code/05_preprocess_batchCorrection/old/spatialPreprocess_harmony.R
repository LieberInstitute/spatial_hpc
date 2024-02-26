# tutorial from: https://edward130603.github.io/BayesSpace/articles/joint_clustering.html

setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
    library("here")
    library("spatialLIBD")
    library("ggplot2")
    library("harmony")
    library("BayesSpace")
    library("schex")
    library("scater")
})

library("patchwork")
library("scran")

# load SPE
load(file = here::here("processed-data", "04_QC", "spe_QC.Rdata"))
dim(spe)
# [1] 30359 137442

spe$discard_auto_id <- spe$low_sum_id | spe$low_detected_id
spe <- spe[, colData(spe)$discard_auto_id == FALSE]
dim(spe)
# [1]  30359 135640

# Pre process
spe <- spatialPreprocess(spe, n.PCs = 50, log.normalize=TRUE, assay.type="logcounts")
# Warning message:
#     In asMethod(object) :
#     sparse->dense coercion: allocating vector of size 2.0 GiB

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

hex <- make_hexbin(spe, nbins = 100, dimension_reduction = "UMAP", use_dims = c(1, 2))

# UMAP plots
pdf(file = here::here("plots", "05_preprocess_batchCorrection", "spatialPreprocess_UMAP.pdf"))
ggplot(
    data.frame(reducedDim(spe, "UMAP")),
    aes(x = UMAP1, y = UMAP2, color = factor(spe$brnum))) +
    geom_point() +
    labs(color = "Subject/Brain") +
    theme_bw()

label_df <- make_hexbin_label(hex, col = "brnum")
plot_hexbin_meta(hex, col = "brnum", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("Brains") + theme(legend.position = "right")

ggplot(
    data.frame(reducedDim(spe, "UMAP")),
    aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))) +
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
# Warning messages:
# 1: Quick-TRANSfer stage steps exceeded maximum (= 6782000) 
# 2: Quick-TRANSfer stage steps exceeded maximum (= 6782000) 
# 3: did not converge in 25 iterations 
# 4: Quick-TRANSfer stage steps exceeded maximum (= 6782000) 
# 5: did not converge in 25 iterations 

spe <- runUMAP(spe, dimred = "HARMONY", name = "UMAP.HARMONY")
# Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
# Also defined by ‘spam’
# Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
# Also defined by ‘spam’

colnames(reducedDim(spe, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")

hex <- make_hexbin(spe, nbins = 100, dimension_reduction = "UMAP.HARMONY", use_dims = c(1, 2))

pdf(file = here::here("plots", "05_preprocess_batchCorrection", "spatialPreprocess_UMAP_harmony.pdf"))
ggplot(
    data.frame(reducedDim(spe, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(spe$brnum))) +
    geom_point() +
    labs(color = "Subject/Brain") +
    theme_bw()

label_df <- make_hexbin_label(hex, col = "brnum")
plot_hexbin_meta(hex, col = "brnum", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("Brains") + theme(legend.position = "right")

ggplot(
    data.frame(reducedDim(spe, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(spe$sample_id))) +
    geom_point() +
    labs(color = "Sample/Capture Area") +
    theme_bw()

label_df <- make_hexbin_label(hex, col = "sample_id")
plot_hexbin_meta(hex, col = "sample_id", action = "majority", xlab = "UMAP1", ylab = "UMAP2") + ggtitle("Capture area") + theme(legend.position = "right")

dev.off()

save(spe, file = here::here("processed-data", "05_preprocess_batchCorrection", "spatialPreprocess_harmony_spe.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()