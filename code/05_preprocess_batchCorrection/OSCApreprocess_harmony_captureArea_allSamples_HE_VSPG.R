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

spe<-spe[, colData(spe)$discard_auto_id == FALSE]
dim(spe)
#[1]  31483 188762
##Compute library factors (OSCA)
set.seed(25323)
spe <- computeLibraryFactors(spe)
summary(sizeFactors(spe))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.02266  0.42134  0.70615  1.00000  1.22312 29.60345
pdf(file = here::here("plots", "05_preprocess_batchCorrection","sizeFactor_histogram.pdf"))
hist(sizeFactors(spe), breaks = 1000,xlim=c(0,10))
dev.off()


message("Running logNormCounts()")
spe <- logNormCounts(spe)

message("Running modelGeneVar()")
## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
dec <- modelGeneVar(spe,
                    block = spe$sample_id#,
                    #BPPARAM = MulticoreParam(4)
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

# message("Running getTopHVGs()")
# top.hvgs <- getTopHVGs(dec, prop = 0.1)
# length(top.hvgs)
# 
# top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
# length(top.hvgs.fdr5)
# 
# top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.00001)
# length(top.hvgs.fdr1)
# 
# save(top.hvgs,
#      top.hvgs.fdr5,
#      top.hvgs.fdr1,
#      file = here::here("processed-data", "05_preprocess_batchCorrection", "top.hvgs_captureArea_allSamples.Rdata")
# )

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

# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-10-09
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
# beachmat               2.16.0    2023-04-25 [2] Bioconductor
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.3.1)
# Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
# BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
# BiocNeighbors          1.18.0    2023-04-25 [2] Bioconductor
# BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
# BiocSingular           1.16.0    2023-04-25 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# bluster                1.10.0    2023-04-25 [2] Bioconductor
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.3.1)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
# DelayedMatrixStats     1.22.6    2023-08-28 [2] Bioconductor
# dplyr                * 1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
# dqrng                  0.3.1     2023-08-30 [2] CRAN (R 4.3.1)
# DropletUtils           1.20.0    2023-04-25 [2] Bioconductor
# edgeR                  3.42.4    2023-05-31 [2] Bioconductor
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
# forcats              * 1.0.0     2023-01-29 [2] CRAN (R 4.3.1)
# fs                     1.6.3     2023-07-20 [2] CRAN (R 4.3.1)
# gargle                 1.5.2     2023-07-20 [2] CRAN (R 4.3.1)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.36.3    2023-09-07 [2] Bioconductor
# GenomeInfoDbData       1.2.10    2023-07-20 [2] Bioconductor
# GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
# ggbeeswarm             0.7.2     2023-04-29 [2] CRAN (R 4.3.1)
# ggplot2              * 3.4.3     2023-08-14 [2] CRAN (R 4.3.1)
# ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.3.1)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.1)
# googledrive            2.1.1     2023-06-11 [2] CRAN (R 4.3.1)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
# HDF5Array              1.28.1    2023-05-01 [2] Bioconductor
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.3.1)
# igraph                 1.5.1     2023-08-10 [2] CRAN (R 4.3.1)
# IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.3.1)
# jaffelab             * 0.99.32   2023-08-18 [1] Github (LieberInstitute/jaffelab@21e6574)
# lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
# limma                  3.56.2    2023-06-04 [2] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
# lubridate            * 1.9.2     2023-02-10 [2] CRAN (R 4.3.1)
# magick                 2.7.5     2023-08-07 [2] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
# MASS                   7.3-60    2023-05-04 [3] CRAN (R 4.3.1)
# Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# metapod                1.8.0     2023-04-25 [2] Bioconductor
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
# nlme                   3.1-163   2023-08-09 [3] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# purrr                * 1.0.2     2023-08-10 [2] CRAN (R 4.3.1)
# R.methodsS3            1.8.2     2022-06-13 [2] CRAN (R 4.3.1)
# R.oo                   1.25.0    2022-06-12 [2] CRAN (R 4.3.1)
# R.utils                2.12.2    2022-11-11 [2] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.3.1)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.3.1)
# Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# readr                * 2.1.4     2023-02-10 [2] CRAN (R 4.3.1)
# rhdf5                  2.44.0    2023-04-25 [2] Bioconductor
# rhdf5filters           1.12.1    2023-04-30 [2] Bioconductor
# Rhdf5lib               1.22.1    2023-09-10 [2] Bioconductor
# rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.3.1)
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.3.1)
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# ScaledMatrix           1.8.1     2023-05-03 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
# scater               * 1.28.0    2023-04-25 [2] Bioconductor
# scran                * 1.28.2    2023-07-23 [2] Bioconductor
# scuttle              * 1.10.2    2023-08-03 [2] Bioconductor
# segmented              1.6-4     2023-04-13 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
# sparseMatrixStats      1.12.2    2023-07-02 [2] Bioconductor
# SpatialExperiment    * 1.10.0    2023-04-25 [2] Bioconductor
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.3.1)
# stringi                1.7.12    2023-01-11 [2] CRAN (R 4.3.1)
# stringr              * 1.5.0     2022-12-02 [2] CRAN (R 4.3.1)
# SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
# tibble               * 3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
# tidyr                * 1.3.0     2023-01-24 [2] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
# tidyverse            * 2.0.0     2023-02-22 [2] CRAN (R 4.3.1)
# timechange             0.2.0     2023-01-11 [2] CRAN (R 4.3.1)
# tzdb                   0.4.0     2023-05-12 [2] CRAN (R 4.3.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
# vctrs                  0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.3.1)
# viridis                0.6.4     2023-07-22 [2] CRAN (R 4.3.1)
# viridisLite            0.4.2     2023-05-02 [2] CRAN (R 4.3.1)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.3.1)
# XVector                0.40.0    2023-04-25 [2] Bioconductor
# zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
# 
# [1] /users/enelson/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
