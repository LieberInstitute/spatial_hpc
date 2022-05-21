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
})

# load SPE
load(file = here::here("processed-data", "pilot_data_checks", "spe.Rdata"))
dim(spe)
# [1] 27633 28871

# Pre process
spePP <- spatialPreprocess(spe, n.PCs = 50)
spePP <- runUMAP(spePP, dimred = "PCA")

# Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
# Also defined by ‘spam’
# Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
# Also defined by ‘spam’
# Warning message:
#   In .check_reddim_names(x, value, withDimnames) :
#   non-NULL 'rownames(value)' should be the same as 'colnames(x)' for
# 'reducedDim<-'. This will be an error in the next release of
# Bioconductor.

colnames(reducedDim(spePP, "UMAP")) <- c("UMAP1", "UMAP2")

# UMAP plots
dir.create(here::here("plots", "pilot_data_checks"), showWarnings = FALSE)
pdf(file = here::here("plots", "pilot_data_checks", "hpc_UMAP_spe.pdf"))
ggplot(
    data.frame(reducedDim(spePP, "UMAP")),
    aes(x = UMAP1, y = UMAP2, color = factor(spePP$sample_id))
) +
    geom_point() +
    labs(color = "Sample/Capture Area") +
    theme_bw()

ggplot(
    data.frame(reducedDim(spePP, "UMAP")),
    aes(x = UMAP1, y = UMAP2, color = factor(spePP$subject))
) +
    geom_point() +
    labs(color = "Subject/Brain") +
    theme_bw()

ggplot(
    data.frame(reducedDim(spePP, "UMAP")),
    aes(x = UMAP1, y = UMAP2, color = factor(spePP$slice))
) +
    geom_point() +
    labs(color = "Tissue Slices") +
    theme_bw()

dev.off()

# batch correction - HARMONY
# install.packages("devtools")
# devtools::install_github("immunogenomics/harmony")

speH <- RunHarmony(spePP, "sample_id", verbose = F)
speH <- runUMAP(speH, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(speH, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")

pdf(file = here::here("plots", "pilot_data_checks", "hpc_UMAP_harmony.pdf"))
ggplot(
    data.frame(reducedDim(speH, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(speH$sample_id))
) +
    geom_point() +
    labs(color = "Sample/Capture Area") +
    theme_bw()

ggplot(
    data.frame(reducedDim(speH, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(speH$subject))
) +
    geom_point() +
    labs(color = "Subject/Brain") +
    theme_bw()

ggplot(
    data.frame(reducedDim(speH, "UMAP.HARMONY")),
    aes(x = UMAP1, y = UMAP2, color = factor(speH$slice))
) +
    geom_point() +
    labs(color = "Tissue Slices") +
    theme_bw()

dev.off()

save(speH, file = here::here("processed-data", "pilot_data_checks", "spe_harmony.Rdata"))
