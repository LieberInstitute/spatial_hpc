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
  library("Polychrome")})


## Load the data
load(file = here::here("processed-data", "04_QC", "spe_QC_allSamples.Rdata"))
dim(spe)
#[1]  31483 191136

###remove discarded spots
spe<-spe[, colData(spe)$discard_auto_id == FALSE]
dim(spe)
#[1]  31483 188762

##Compute library factors
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
dec <- modelGeneVar(spe,
                    block = spe$sample_id#,
                    #BPPARAM = MulticoreParam(4)
)


 message("Running getTopHVGs()")
 top.hvgs <- getTopHVGs(dec, prop = 0.1)
 length(top.hvgs)
save(top.hvgs,
     file = here::here("processed-data", "05_preprocess_batchCorrection", "top.hvgs_captureArea_allSamples.Rdata")
)

####dimred
 message("Running runPCA()")
 set.seed(20220201)
 Sys.time()
 spe <- runPCA(spe, subset_row = top.hvgs, ncomponents = 50)
 Sys.time()

 message("Running runUMAP()")
 Sys.time()
 set.seed(20220208)
 spe <- runUMAP(spe, dimred = "PCA")
 colnames(reducedDim(spe, "UMAP")) <- c("UMAP1", "UMAP2")
 Sys.time()

####make fig S10 (batch effect UMAP)
spe$VSPG<-ifelse(spe$brnum %in% levels(spe$brnum)[11:12],T,F)

pdf(file=here::here('plots','figures','supp_figures','fig_s10','fig_s10.pdf'))
plotUMAP(spe,colour_by='VSPG',point_size=0.5,theme_size=20)
dev.off()

save(spe,file=here::here('processed-data','05_preprocess_batchCorrection','spe_norm_final.rda'))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

#[1] "Reproducibility information:"
#[1] "2024-02-25 22:53:12 EST"
#     user    system   elapsed
#  745.203    19.195 27790.656
#- Session info ---------------------------------------------------------------------------------------
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-02-25
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
#
#- Packages -------------------------------------------------------------------------------------------
# package                * version   date (UTC) lib source
# abind                    1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
# AnnotationDbi            1.62.2    2023-07-02 [2] Bioconductor
# AnnotationHub            3.8.0     2023-04-25 [2] Bioconductor
# attempt                  0.3.1     2020-05-03 [2] CRAN (R 4.3.1)
# beachmat                 2.16.0    2023-04-25 [2] Bioconductor
# beeswarm                 0.4.0     2021-06-01 [2] CRAN (R 4.3.1)
# benchmarkme              1.0.8     2022-06-12 [2] CRAN (R 4.3.1)
# benchmarkmeData          1.0.4     2020-04-23 [2] CRAN (R 4.3.1)
# Biobase                * 2.60.0    2023-04-25 [2] Bioconductor
# BiocFileCache            2.8.0     2023-04-25 [2] Bioconductor
# BiocGenerics           * 0.46.0    2023-04-25 [2] Bioconductor
# BiocIO                   1.10.0    2023-04-25 [2] Bioconductor
# BiocManager              1.30.22   2023-08-08 [2] CRAN (R 4.3.1)
# BiocNeighbors            1.18.0    2023-04-25 [2] Bioconductor
# BiocParallel           * 1.34.2    2023-05-22 [2] Bioconductor
# BiocSingular             1.16.0    2023-04-25 [2] Bioconductor
# BiocVersion              3.17.1    2022-11-04 [2] Bioconductor
# Biostrings               2.68.1    2023-05-16 [2] Bioconductor
# bit                      4.0.5     2022-11-15 [2] CRAN (R 4.3.1)
# bit64                    4.0.5     2020-08-30 [2] CRAN (R 4.3.1)
# bitops                   1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# blob                     1.2.4     2023-03-17 [2] CRAN (R 4.3.1)
# bluster                  1.10.0    2023-04-25 [2] Bioconductor
# bslib                    0.5.1     2023-08-11 [2] CRAN (R 4.3.1)
# cachem                   1.0.8     2023-05-01 [2] CRAN (R 4.3.1)
# cli                      3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# cluster                  2.1.4     2022-08-22 [3] CRAN (R 4.3.1)
# codetools                0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
# colorspace               2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# config                   0.3.2     2023-08-30 [2] CRAN (R 4.3.1)
# cowplot                  1.1.1     2020-12-30 [2] CRAN (R 4.3.1)
# crayon                   1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# curl                     5.0.2     2023-08-14 [2] CRAN (R 4.3.1)
# data.table               1.14.8    2023-02-17 [2] CRAN (R 4.3.1)
# DBI                      1.1.3     2022-06-18 [2] CRAN (R 4.3.1)
# dbplyr                   2.3.3     2023-07-07 [2] CRAN (R 4.3.1)
# DelayedArray             0.26.7    2023-07-28 [2] Bioconductor
# DelayedMatrixStats       1.22.6    2023-08-28 [2] Bioconductor
# digest                   0.6.33    2023-07-07 [2] CRAN (R 4.3.1)
# doParallel               1.0.17    2022-02-07 [2] CRAN (R 4.3.1)
# dotCall64                1.0-2     2022-10-03 [2] CRAN (R 4.3.1)
# dplyr                  * 1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
# dqrng                    0.3.1     2023-08-30 [2] CRAN (R 4.3.1)
# DropletUtils             1.20.0    2023-04-25 [2] Bioconductor
# DT                       0.29      2023-08-29 [2] CRAN (R 4.3.1)
# edgeR                    3.42.4    2023-05-31 [2] Bioconductor
# ellipsis                 0.3.2     2021-04-29 [2] CRAN (R 4.3.1)
# ExperimentHub            2.8.1     2023-07-12 [2] Bioconductor
# fansi                    1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
# farver                   2.1.1     2022-07-06 [2] CRAN (R 4.3.1)
# fastmap                  1.1.1     2023-02-24 [2] CRAN (R 4.3.1)
# fields                   15.2      2023-08-17 [2] CRAN (R 4.3.1)
# filelock                 1.0.2     2018-10-05 [2] CRAN (R 4.3.1)
# forcats                * 1.0.0     2023-01-29 [2] CRAN (R 4.3.1)
# foreach                  1.5.2     2022-02-02 [2] CRAN (R 4.3.1)
# fs                       1.6.3     2023-07-20 [2] CRAN (R 4.3.1)
# gargle                   1.5.2     2023-07-20 [2] CRAN (R 4.3.1)
# generics                 0.1.3     2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb           * 1.36.3    2023-09-07 [2] Bioconductor
# GenomeInfoDbData         1.2.10    2023-07-20 [2] Bioconductor
# GenomicAlignments        1.36.0    2023-04-25 [2] Bioconductor
# GenomicRanges          * 1.52.0    2023-04-25 [2] Bioconductor
# ggbeeswarm               0.7.2     2023-04-29 [2] CRAN (R 4.3.1)
# ggplot2                * 3.4.3     2023-08-14 [2] CRAN (R 4.3.1)
# ggrepel                * 0.9.3     2023-02-03 [2] CRAN (R 4.3.1)
# glue                     1.6.2     2022-02-24 [2] CRAN (R 4.3.1)
# golem                    0.4.1     2023-06-05 [2] CRAN (R 4.3.1)
# googledrive              2.1.1     2023-06-11 [2] CRAN (R 4.3.1)
# gridExtra              * 2.3       2017-09-09 [2] CRAN (R 4.3.1)
# gtable                   0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
# harmony                * 1.0.1     2023-09-20 [2] CRAN (R 4.3.1)
# HDF5Array                1.28.1    2023-05-01 [2] Bioconductor
# here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# hms                      1.1.3     2023-03-21 [2] CRAN (R 4.3.1)
# htmltools                0.5.6     2023-08-10 [2] CRAN (R 4.3.1)
# htmlwidgets              1.6.2     2023-03-17 [2] CRAN (R 4.3.1)
# httpuv                   1.6.11    2023-05-11 [2] CRAN (R 4.3.1)
# httr                     1.4.7     2023-08-15 [2] CRAN (R 4.3.1)
# igraph                   1.5.1     2023-08-10 [2] CRAN (R 4.3.1)
# interactiveDisplayBase   1.38.0    2023-04-25 [2] Bioconductor
# IRanges                * 2.34.1    2023-06-22 [2] Bioconductor
# irlba                    2.3.5.1   2022-10-03 [2] CRAN (R 4.3.1)
# iterators                1.0.14    2022-02-05 [2] CRAN (R 4.3.1)
# jaffelab               * 0.99.32   2023-08-18 [1] Github (LieberInstitute/jaffelab@21e6574)
# jquerylib                0.1.4     2021-04-26 [2] CRAN (R 4.3.1)
# jsonlite                 1.8.7     2023-06-29 [2] CRAN (R 4.3.1)
# KEGGREST                 1.40.0    2023-04-25 [2] Bioconductor
# labeling                 0.4.3     2023-08-29 [2] CRAN (R 4.3.1)
# later                    1.3.1     2023-05-02 [2] CRAN (R 4.3.1)
# lattice                  0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# lazyeval                 0.2.2     2019-03-15 [2] CRAN (R 4.3.1)
# lifecycle                1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
# limma                    3.56.2    2023-06-04 [2] Bioconductor
# locfit                   1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
# lubridate              * 1.9.2     2023-02-10 [2] CRAN (R 4.3.1)
# magick                   2.7.5     2023-08-07 [2] CRAN (R 4.3.1)
# magrittr                 2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
# maps                     3.4.1     2022-10-30 [2] CRAN (R 4.3.1)
# MASS                     7.3-60    2023-05-04 [3] CRAN (R 4.3.1)
# Matrix                 * 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics         * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats            * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# memoise                  2.0.1     2021-11-26 [2] CRAN (R 4.3.1)
# metapod                  1.8.0     2023-04-25 [2] Bioconductor
# mime                     0.12      2021-09-28 [2] CRAN (R 4.3.1)
# munsell                  0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
# nlme                     3.1-163   2023-08-09 [3] CRAN (R 4.3.1)
# paletteer                1.5.0     2022-10-19 [2] CRAN (R 4.3.1)
# PCAtools               * 2.12.0    2023-04-25 [1] Bioconductor
# pillar                   1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# plotly                   4.10.2    2023-06-03 [2] CRAN (R 4.3.1)
# plyr                     1.8.8     2022-11-11 [2] CRAN (R 4.3.1)
# png                      0.1-8     2022-11-29 [2] CRAN (R 4.3.1)
# Polychrome             * 1.5.1     2022-05-03 [1] CRAN (R 4.3.1)
# promises                 1.2.1     2023-08-10 [2] CRAN (R 4.3.1)
# purrr                  * 1.0.2     2023-08-10 [2] CRAN (R 4.3.1)
# R.methodsS3              1.8.2     2022-06-13 [2] CRAN (R 4.3.1)
# R.oo                     1.25.0    2022-06-12 [2] CRAN (R 4.3.1)
# R.utils                  2.12.2    2022-11-11 [2] CRAN (R 4.3.1)
# R6                       2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# rafalib                * 1.0.0     2015-08-09 [1] CRAN (R 4.3.1)
# rappdirs                 0.3.3     2021-01-31 [2] CRAN (R 4.3.1)
# RColorBrewer             1.1-3     2022-04-03 [2] CRAN (R 4.3.1)
# Rcpp                   * 1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RcppAnnoy                0.0.21    2023-07-02 [2] CRAN (R 4.3.1)
# RCurl                    1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# readr                  * 2.1.4     2023-02-10 [2] CRAN (R 4.3.1)
# rematch2                 2.1.2     2020-05-01 [2] CRAN (R 4.3.1)
# reshape2                 1.4.4     2020-04-09 [2] CRAN (R 4.3.1)
# restfulr                 0.0.15    2022-06-16 [2] CRAN (R 4.3.1)
# rhdf5                    2.44.0    2023-04-25 [2] Bioconductor
# rhdf5filters             1.12.1    2023-04-30 [2] Bioconductor
# Rhdf5lib                 1.22.1    2023-09-10 [2] Bioconductor
# rjson                    0.2.21    2022-01-09 [2] CRAN (R 4.3.1)
# rlang                    1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot                2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# Rsamtools                2.16.0    2023-04-25 [2] Bioconductor
# RSQLite                  2.3.1     2023-04-03 [2] CRAN (R 4.3.1)
# rsvd                     1.0.5     2021-04-16 [2] CRAN (R 4.3.1)
# rtracklayer              1.60.1    2023-08-15 [2] Bioconductor
# S4Arrays                 1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors              * 0.38.1    2023-05-02 [2] Bioconductor
# sass                     0.4.7     2023-07-15 [2] CRAN (R 4.3.1)
# ScaledMatrix             1.8.1     2023-05-03 [2] Bioconductor
# scales                   1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
# scater                 * 1.28.0    2023-04-25 [2] Bioconductor
# scatterplot3d            0.3-44    2023-05-05 [1] CRAN (R 4.3.1)
# scran                  * 1.28.2    2023-07-23 [2] Bioconductor
# scuttle                * 1.10.2    2023-08-03 [2] Bioconductor
# segmented                1.6-4     2023-04-13 [1] CRAN (R 4.3.1)
# sessioninfo            * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# shiny                    1.7.5     2023-08-12 [2] CRAN (R 4.3.1)
# shinyWidgets             0.8.0     2023-08-30 [2] CRAN (R 4.3.1)
# SingleCellExperiment   * 1.22.0    2023-04-25 [2] Bioconductor
# spam                     2.9-1     2022-08-07 [2] CRAN (R 4.3.1)
# sparseMatrixStats        1.12.2    2023-07-02 [2] Bioconductor
# SpatialExperiment      * 1.10.0    2023-04-25 [2] Bioconductor
# spatialLIBD            * 1.12.0    2023-04-27 [2] Bioconductor
# statmod                  1.5.0     2023-01-06 [2] CRAN (R 4.3.1)
# stringi                  1.7.12    2023-01-11 [2] CRAN (R 4.3.1)
# stringr                * 1.5.0     2022-12-02 [2] CRAN (R 4.3.1)
# SummarizedExperiment   * 1.30.2    2023-06-06 [2] Bioconductor
# tibble                 * 3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
# tidyr                  * 1.3.0     2023-01-24 [2] CRAN (R 4.3.1)
# tidyselect               1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
# tidyverse              * 2.0.0     2023-02-22 [2] CRAN (R 4.3.1)
# timechange               0.2.0     2023-01-11 [2] CRAN (R 4.3.1)
# tzdb                     0.4.0     2023-05-12 [2] CRAN (R 4.3.1)
# utf8                     1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
# uwot                   * 0.1.16    2023-06-29 [2] CRAN (R 4.3.1)
# vctrs                    0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
# vipor                    0.4.5     2017-03-22 [2] CRAN (R 4.3.1)
# viridis                  0.6.4     2023-07-22 [2] CRAN (R 4.3.1)
# viridisLite              0.4.2     2023-05-02 [2] CRAN (R 4.3.1)
# withr                    2.5.0     2022-03-03 [2] CRAN (R 4.3.1)
# XML                      3.99-0.14 2023-03-19 [2] CRAN (R 4.3.1)
# xtable                   1.8-4     2019-04-21 [2] CRAN (R 4.3.1)
# XVector                  0.40.0    2023-04-25 [2] Bioconductor
# yaml                     2.3.7     2023-01-23 [2] CRAN (R 4.3.1)
# zlibbioc                 1.46.0    2023-04-25 [2] Bioconductor
#
# [1] /users/enelson/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
#
#------------------------------------------------------------------------------------------------------


