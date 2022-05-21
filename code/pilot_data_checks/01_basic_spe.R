
# cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/
suppressPackageStartupMessages(library("here"))
# remotes::install_github("drighelli/SpatialExperiment")
# remotes::install_github("LieberInstitute/spatialLIBD")
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("lobstr"))
suppressPackageStartupMessages(library("sessioninfo"))

## Define some info for the samples
sample_info <- data.frame(
    subject = c("Br6423-O", "Br6423-O", "Br6423-O", "Br6423-O", "Br6432-R", "Br6432-R", "Br2743-Y", "Br2743-Y"),
    slice = c("Br6423-O_1", "Br6423-O_1", "Br6423-O_2", "Br6423-O_2", "Br6432-R_1", "Br6432-R_1", "Br2743-Y_1", "Br2743-Y_1"),
    sample_id = c(
        "V10B01-085_A1",
        "V10B01-085_B1",
        "V10B01-085_C1",
        "V10B01-085_D1",
        "V10B01-086_A1",
        "V10B01-086_B1",
        "V10B01-086_C1",
        "V10B01-086_D1"
    )
)
sample_info$sample_path <-
    file.path(
        here::here("processed-data", "spaceranger_novaseq"),
        sample_info$sample_id,
        "outs"
    )
stopifnot(all(file.exists(sample_info$sample_path)))

## Define the donor info using information from
## https://github.com/LieberInstitute/spatial_hpc/blob/main/raw-data/sample_info/Visium_HPC_Round1%2B2_110321_Master_SCP.xlsx
donor_info <- data.frame(
    subject = c("Br6423-O", "Br6432-R", "Br2743-Y"),
    age = c(51.73, 48.88, 61.54),
    sex = c("M", "M", "M"),
    # race = "EA/CAUC",
    diagnosis = "Control"
    # rin = c(),
    # pmi = c(),
)

## Combine sample info with the donor info
sample_info <- merge(sample_info, donor_info)

## Build basic SPE
Sys.time()
spe <- read10xVisiumWrapper(
    sample_info$sample_path,
    sample_info$sample_id,
    type = "sparse",
    data = "raw",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE
)
Sys.time()

# [1] "2022-03-04 15:19:02 EST"
# 2022-03-04 15:19:05 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
# 2022-03-04 15:22:15 read10xVisiumAnalysis: reading analysis output from SpaceRanger
# 2022-03-04 15:22:28 add10xVisiumAnalysis: adding analysis output from SpaceRanger
# 2022-03-04 15:22:33 rtracklayer::import: reading the reference GTF file
# 2022-03-04 15:24:35 adding gene information to the SPE object
# 2022-03-04 15:24:35 adding information used by spatialLIBD
# [1] "2022-03-04 15:24:50 EST"

## Add the study design info
add_design <- function(spe) {
    new_col <- merge(colData(spe), sample_info)
    ## Fix order
    new_col <- new_col[match(spe$key, new_col$key), ]
    stopifnot(identical(new_col$key, spe$key))
    rownames(new_col) <- rownames(colData(spe))
    colData(spe) <-
        new_col[, -which(colnames(new_col) == "sample_path")]
    return(spe)
}
spe <- add_design(spe)

## Read in cell counts and segmentation results
segmentations_list <-
    lapply(sample_info$sample_id, function(sampleid) {
        file <-
            here(
                "processed-data",
                "spaceranger_novaseq",
                sampleid,
                "outs",
                "spatial",
                "tissue_spot_counts.csv"
            )
        if (!file.exists(file)) {
            return(NULL)
        }
        x <- read.csv(file)
        x$key <- paste0(x$barcode, "_", sampleid)
        return(x)
    })

## Merge them (once the these files are done, this could be replaced by an rbind)
segmentations <-
    Reduce(function(...) {
        merge(..., all = TRUE)
    }, segmentations_list[lengths(segmentations_list) > 0])

## Add the information
segmentation_match <- match(spe$key, segmentations$key)
segmentation_info <-
    segmentations[segmentation_match, -which(
        colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
    )]
colData(spe) <- cbind(colData(spe), segmentation_info)

## basic SPE
spe_basic <- spe
# dir.create(here::here("processed-data", "pilot_data_checks"), showWarnings = FALSE)
save(spe_basic, file = here::here("processed-data", "pilot_data_checks", "spe_basic.Rdata"))

## Size in Gb
lobstr::obj_size(spe_basic) / 1024^3
# 1.344729 B
dim(spe_basic)
# [1] 36601 39936

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


# Reproducibility information
# [1] "2022-03-04 15:45:31 EST"
# user   system  elapsed
# 853.217   75.442 2830.776
#
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-03-04
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version  date (UTC) lib source
# AnnotationDbi            1.56.2   2021-11-09 [2] Bioconductor
# AnnotationHub            3.2.2    2022-03-01 [2] Bioconductor
# assertthat               0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# attempt                  0.3.1    2020-05-03 [1] CRAN (R 4.1.1)
# beachmat                 2.10.0   2021-10-26 [2] Bioconductor
# beeswarm                 0.4.0    2021-06-01 [1] CRAN (R 4.1.1)
# benchmarkme              1.0.7    2021-03-21 [1] CRAN (R 4.1.1)
# benchmarkmeData          1.0.4    2020-04-23 [1] CRAN (R 4.1.1)
# Biobase                * 2.54.0   2021-10-26 [2] Bioconductor
# BiocFileCache            2.2.1    2022-01-23 [2] Bioconductor
# BiocGenerics           * 0.40.0   2021-10-26 [1] Bioconductor
# BiocIO                   1.4.0    2021-10-26 [2] Bioconductor
# BiocManager              1.30.16  2021-06-15 [1] CRAN (R 4.1.2)
# BiocNeighbors            1.12.0   2021-10-26 [1] Bioconductor
# BiocParallel             1.28.3   2021-12-09 [1] Bioconductor
# BiocSingular             1.10.0   2021-10-26 [1] Bioconductor
# BiocVersion              3.14.0   2021-05-19 [2] Bioconductor
# Biostrings               2.62.0   2021-10-26 [2] Bioconductor
# bit                      4.0.4    2020-08-04 [2] CRAN (R 4.1.0)
# bit64                    4.0.5    2020-08-30 [2] CRAN (R 4.1.0)
# bitops                   1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# blob                     1.2.2    2021-07-23 [2] CRAN (R 4.1.0)
# brio                     1.1.3    2021-11-30 [2] CRAN (R 4.1.2)
# bslib                    0.3.1    2021-10-06 [2] CRAN (R 4.1.2)
# cachem                   1.0.6    2021-08-19 [2] CRAN (R 4.1.2)
# callr                    3.7.0    2021-04-20 [2] CRAN (R 4.1.0)
# cli                      3.2.0    2022-02-14 [2] CRAN (R 4.1.2)
# codetools                0.2-18   2020-11-04 [3] CRAN (R 4.1.2)
# colorout               * 1.2-2    2022-03-04 [1] Github (jalvesaq/colorout@79931fd)
# colorspace               2.0-3    2022-02-21 [2] CRAN (R 4.1.2)
# config                   0.3.1    2020-12-17 [1] CRAN (R 4.1.1)
# cowplot                  1.1.1    2020-12-30 [1] CRAN (R 4.1.1)
# crayon                   1.5.0    2022-02-14 [2] CRAN (R 4.1.2)
# curl                     4.3.2    2021-06-23 [2] CRAN (R 4.1.0)
# data.table               1.14.2   2021-09-27 [1] CRAN (R 4.1.2)
# DBI                      1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
# dbplyr                   2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
# DelayedArray             0.20.0   2021-10-26 [1] Bioconductor
# DelayedMatrixStats       1.16.0   2021-10-26 [1] Bioconductor
# desc                     1.4.0    2021-09-28 [2] CRAN (R 4.1.2)
# digest                   0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
# doParallel               1.0.17   2022-02-07 [2] CRAN (R 4.1.2)
# dotCall64                1.0-1    2021-02-11 [2] CRAN (R 4.1.0)
# dplyr                    1.0.8    2022-02-08 [2] CRAN (R 4.1.2)
# dqrng                    0.3.0    2021-05-01 [1] CRAN (R 4.1.1)
# DropletUtils             1.14.2   2022-01-09 [1] Bioconductor
# DT                       0.21     2022-02-26 [2] CRAN (R 4.1.2)
# edgeR                    3.36.0   2021-10-26 [1] Bioconductor
# ellipsis                 0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# ExperimentHub            2.2.1    2022-01-23 [2] Bioconductor
# fansi                    1.0.2    2022-01-14 [2] CRAN (R 4.1.2)
# fastmap                  1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
# fields                   13.3     2021-10-30 [2] CRAN (R 4.1.2)
# filelock                 1.0.2    2018-10-05 [2] CRAN (R 4.1.0)
# foreach                  1.5.2    2022-02-02 [2] CRAN (R 4.1.2)
# fs                       1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
# generics                 0.1.2    2022-01-31 [2] CRAN (R 4.1.2)
# GenomeInfoDb           * 1.30.1   2022-01-30 [1] Bioconductor
# GenomeInfoDbData         1.2.7    2021-11-01 [2] Bioconductor
# GenomicAlignments        1.30.0   2021-10-26 [2] Bioconductor
# GenomicRanges          * 1.46.1   2021-11-18 [2] Bioconductor
# ggbeeswarm               0.6.0    2017-08-07 [1] CRAN (R 4.1.1)
# ggplot2                  3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
# ggrepel                  0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
# glue                     1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
# golem                    0.3.2    2022-03-04 [1] CRAN (R 4.1.2)
# gridExtra                2.3      2017-09-09 [2] CRAN (R 4.1.0)
# gtable                   0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# HDF5Array                1.22.1   2021-11-14 [2] Bioconductor
# here                   * 1.0.1    2020-12-13 [1] CRAN (R 4.1.1)
# htmltools                0.5.2    2021-08-25 [2] CRAN (R 4.1.2)
# htmlwidgets              1.5.4    2021-09-08 [2] CRAN (R 4.1.2)
# httpuv                   1.6.5    2022-01-05 [2] CRAN (R 4.1.2)
# httr                     1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
# interactiveDisplayBase   1.32.0   2021-10-26 [2] Bioconductor
# IRanges                * 2.28.0   2021-10-26 [1] Bioconductor
# irlba                    2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
# iterators                1.0.14   2022-02-05 [2] CRAN (R 4.1.2)
# jquerylib                0.1.4    2021-04-26 [2] CRAN (R 4.1.0)
# jsonlite                 1.8.0    2022-02-22 [2] CRAN (R 4.1.2)
# KEGGREST                 1.34.0   2021-10-26 [2] Bioconductor
# knitr                    1.37     2021-12-16 [2] CRAN (R 4.1.2)
# later                    1.3.0    2021-08-18 [2] CRAN (R 4.1.2)
# lattice                  0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lazyeval                 0.2.2    2019-03-15 [2] CRAN (R 4.1.0)
# lifecycle                1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
# limma                    3.50.1   2022-02-17 [2] Bioconductor
# lobstr                 * 1.1.1    2019-07-02 [2] CRAN (R 4.1.0)
# locfit                   1.5-9.5  2022-03-03 [2] CRAN (R 4.1.2)
# magick                   2.7.3    2021-08-18 [2] CRAN (R 4.1.2)
# magrittr                 2.0.2    2022-01-26 [2] CRAN (R 4.1.2)
# maps                     3.4.0    2021-09-25 [2] CRAN (R 4.1.2)
# Matrix                   1.4-0    2021-12-08 [3] CRAN (R 4.1.2)
# MatrixGenerics         * 1.6.0    2021-10-26 [1] Bioconductor
# matrixStats            * 0.61.0   2021-09-17 [1] CRAN (R 4.1.2)
# memoise                  2.0.1    2021-11-26 [2] CRAN (R 4.1.2)
# mime                     0.12     2021-09-28 [2] CRAN (R 4.1.2)
# munsell                  0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# pillar                   1.7.0    2022-02-01 [2] CRAN (R 4.1.2)
# pkgbuild                 1.3.1    2021-12-20 [2] CRAN (R 4.1.2)
# pkgconfig                2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# pkgload                  1.2.4    2021-11-30 [2] CRAN (R 4.1.2)
# plotly                   4.10.0   2021-10-09 [2] CRAN (R 4.1.2)
# png                      0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
# Polychrome               1.3.1    2021-07-16 [1] CRAN (R 4.1.1)
# prettyunits              1.1.1    2020-01-24 [2] CRAN (R 4.1.0)
# processx                 3.5.2    2021-04-30 [2] CRAN (R 4.1.0)
# promises                 1.2.0.1  2021-02-11 [2] CRAN (R 4.1.0)
# ps                       1.6.0    2021-02-28 [2] CRAN (R 4.1.0)
# purrr                    0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R.methodsS3              1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
# R.oo                     1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
# R.utils                  2.11.0   2021-09-26 [2] CRAN (R 4.1.2)
# R6                       2.5.1    2021-08-19 [1] CRAN (R 4.1.1)
# rappdirs                 0.3.3    2021-01-31 [2] CRAN (R 4.1.0)
# RColorBrewer             1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                     1.0.8    2022-01-13 [2] CRAN (R 4.1.2)
# RCurl                    1.98-1.6 2022-02-08 [1] CRAN (R 4.1.2)
# remotes                  2.4.2    2021-11-30 [1] CRAN (R 4.1.2)
# restfulr                 0.0.13   2017-08-06 [2] CRAN (R 4.1.0)
# rhdf5                    2.38.0   2021-10-26 [2] Bioconductor
# rhdf5filters             1.6.0    2021-10-26 [2] Bioconductor
# Rhdf5lib                 1.16.0   2021-10-26 [2] Bioconductor
# rjson                    0.2.21   2022-01-09 [2] CRAN (R 4.1.2)
# rlang                    1.0.2    2022-03-04 [1] CRAN (R 4.1.2)
# roxygen2                 7.1.2    2021-09-08 [2] CRAN (R 4.1.2)
# rprojroot                2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# Rsamtools                2.10.0   2021-10-26 [2] Bioconductor
# RSQLite                  2.2.10   2022-02-17 [2] CRAN (R 4.1.2)
# rstudioapi               0.13     2020-11-12 [2] CRAN (R 4.1.0)
# rsvd                     1.0.5    2021-04-16 [1] CRAN (R 4.1.1)
# rtracklayer            * 1.54.0   2021-10-26 [2] Bioconductor
# S4Vectors              * 0.32.3   2021-11-21 [1] Bioconductor
# sass                     0.4.0    2021-05-12 [2] CRAN (R 4.1.0)
# ScaledMatrix             1.2.0    2021-10-26 [1] Bioconductor
# scales                   1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# scater                   1.22.0   2021-10-26 [1] Bioconductor
# scatterplot3d            0.3-41   2018-03-14 [1] CRAN (R 4.1.1)
# scuttle                  1.4.0    2021-10-26 [1] Bioconductor
# sessioninfo            * 1.2.2    2021-12-06 [1] CRAN (R 4.1.2)
# shiny                    1.7.1    2021-10-02 [2] CRAN (R 4.1.2)
# shinyWidgets             0.6.4    2022-02-06 [1] CRAN (R 4.1.2)
# SingleCellExperiment   * 1.16.0   2021-10-26 [1] Bioconductor
# spam                     2.8-0    2022-01-06 [2] CRAN (R 4.1.2)
# sparseMatrixStats        1.6.0    2021-10-26 [1] Bioconductor
# SpatialExperiment      * 1.5.3    2022-03-04 [1] Github (drighelli/SpatialExperiment@caa4209)
# spatialLIBD            * 1.7.13   2022-03-04 [1] Github (LieberInstitute/spatialLIBD@c6a0e18)
# stringi                  1.7.6    2021-11-29 [1] CRAN (R 4.1.2)
# stringr                  1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment   * 1.24.0   2021-10-26 [1] Bioconductor
# testthat                 3.1.2    2022-01-20 [2] CRAN (R 4.1.2)
# tibble                   3.1.6    2021-11-07 [1] CRAN (R 4.1.2)
# tidyr                    1.2.0    2022-02-01 [2] CRAN (R 4.1.2)
# tidyselect               1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
# usethis                  2.1.5    2021-12-09 [2] CRAN (R 4.1.2)
# utf8                     1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                    0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# vipor                    0.4.5    2017-03-22 [1] CRAN (R 4.1.1)
# viridis                  0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
# viridisLite              0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
# withr                    2.5.0    2022-03-03 [2] CRAN (R 4.1.2)
# xfun                     0.30     2022-03-02 [2] CRAN (R 4.1.2)
# XML                      3.99-0.9 2022-02-24 [2] CRAN (R 4.1.2)
# xml2                     1.3.3    2021-11-30 [2] CRAN (R 4.1.2)
# xtable                   1.8-4    2019-04-21 [2] CRAN (R 4.1.0)
# XVector                  0.34.0   2021-10-26 [2] Bioconductor
# yaml                     2.3.5    2022-02-21 [2] CRAN (R 4.1.2)
# zlibbioc                 1.40.0   2021-10-26 [2] Bioconductor
#
# [1] /users/mtippani/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
