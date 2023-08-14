#srun --partition=interactive --mem-per-cpu=12G --cpus-per-task=4 --pty --x11 bash
setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
library("here")
library("SpatialExperiment")
library("spatialLIBD")
library("rtracklayer")
library("lobstr")
library("sessioninfo")
})

## Define some info for the samples
load(here::here("code", "REDCap", "REDCap_HPC.rda"))

sample_info <- data.frame(dateImg = as.Date(REDCap_HPC$date)) 
sample_info$experimenterImg <- as.factor(REDCap_HPC$experimenter_img)
sample_info$slide <- as.factor(REDCap_HPC$slide)
sample_info$array <- as.factor(REDCap_HPC$array)
sample_info$brnum <- as.factor(sapply(strsplit(REDCap_HPC$sample, "-"), `[`, 1))
sample_info$position <- as.factor(REDCap_HPC$adjacent)
sample_info$seqNum <- as.factor(REDCap_HPC$sample_number)
sample_info$experimenterSeq <- as.factor(REDCap_HPC$experimenter_seq)
sample_info$sample_id <- paste(sample_info$slide, sample_info$array, sep = "_")

sample_info$sample_path[sample_info$dateImg <= "2021-10-11"] <- file.path(here::here("processed-data", "01_spaceranger", "spaceranger_novaseq"), sample_info$sample_id, "outs")
sample_info$sample_path[sample_info$dateImg > "2021-10-11"] <- file.path(here::here("processed-data", "01_spaceranger", "spaceranger_2022-04-12_SPag033122"), sample_info$sample_id[sample_info$dateImg > "2021-10-11"], "outs")
stopifnot(all(file.exists(sample_info$sample_path)))

## Define the donor info using information from
donor_info <- read.csv(file.path(here::here("raw-data", "sample_info", "demographicInfo_Geo.csv")), header = TRUE, stringsAsFactors = FALSE)
donor_info <- donor_info[-1]

## Combine sample info with the donor info
sample_info <- merge(sample_info, donor_info)

###We gotta make the round 9 info object too!
sample_info9<-data.frame(
  dateImg=rep(as.Date('2023-01-31'),4),
  experimenterImg=factor(rep("Stephanie Page",4)),
  slide=factor(rep('V12F14-051',4)),
  array=factor(c('A1','B1','C1','D1')),
  brnum=factor(rep('Br2720',4)),
  position=factor(c('BL','BR','TL','TR')),
  seqNum=factor(c("33-10v", "34v_scp", "35v_scp","36v_scp")),
  experimenterSeq=factor(rep("Stephanie Page",4)),
  sample_id = c(
    "Br2720_A1",
    "Br2720_B1",
    "Br2720_C1",
    "Br2720_D1"
  )
)

sample_info9$sample_path <-
  file.path(
    here::here("processed-data", "01_spaceranger", "spaceranger_2023-01-31_round9"),
    sample_info9$sample_id,
    "outs"
  )
stopifnot(all(file.exists(sample_info9$sample_path)))

donor_info9 <- data.frame(
  age = c(48.2, 48.2, 48.2, 48.2),
  sex = c("F", "F", "F", "F"),
  race = c("CAUC", "CAUC", "CAUC", "CAUC"),
  dx = c("Control", "Control", "Control", "Control"),
  pmi = c(25.5, 25.5, 25.5, 25.5)
)

sample_info9 <- cbind(sample_info9, donor_info9)

##match up the colnames
sample_info9<-sample_info9[,match(colnames(sample_info),colnames(sample_info9))]

##rbind the sample_infos
sample_info<-rbind(sample_info,sample_info9)

###We gotta make the round 10 info object too!
sample_info10<-data.frame(
  dateImg=rep(as.Date('2023-06-26'),8),
  experimenterImg=factor(rep("Stephanie Page",8)),
  slide=factor(c(rep('V12D07-332',4),rep('V12D07-335',4))),
  array=factor(c('A1','B1','C1','D1',
                 'A1','B1','C1','D1')),
  brnum=factor(c(rep('Br3942',4),rep('Br8325',4))),
  position=factor(c('TL','BR','BL','TR',
                    'BR','TL','BL','TR')),
  seqNum=factor(paste0(c(51:58),'v_hpc')),
  experimenterSeq=factor(rep("Stephanie Page",8)),
  sample_id = c(
    paste0('V12D07-332_',c('A1','B1','C1','D1'),'_Br3942'),
    paste0('V12D07-335_',c('A1','B1','C1','D1'),'_Br8325')
  )
)

sample_info10$sample_path <-
  file.path(
    here::here("processed-data", "01_spaceranger", "spaceranger_VSPG"),
    sample_info10$sample_id,
    "outs"
  )
stopifnot(all(file.exists(sample_info10$sample_path)))

donor_info10 <- donor_info[match(sample_info10$brnum,donor_info$brnum),]

sample_info10 <- cbind(sample_info10, donor_info10)

##match up the colnames
sample_info10<-sample_info10[,match(colnames(sample_info),colnames(sample_info10))]

##rbind the sample_infos
sample_info<-rbind(sample_info,sample_info10)

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


# dir.create(here::here("processed-data", "pilot_data_checks"), showWarnings = FALSE)
save(spe, file = here::here("processed-data", "02_build_spe", "spe_raw_allSamples.Rdata"))

## Size in Gb
lobstr::obj_size(spe)
# 8.03 GB
dim(spe)
#[1]  36601 219648

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2023-08-11 12:54:36 EDT"
# user   system  elapsed
# 1001.079   30.946 5540.024
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-08-11
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────
# package                * version   date (UTC) lib source
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
# BiocManager              1.30.21.1 2023-07-18 [2] CRAN (R 4.3.1)
# BiocNeighbors            1.18.0    2023-04-25 [2] Bioconductor
# BiocParallel             1.34.2    2023-05-22 [2] Bioconductor
# BiocSingular             1.16.0    2023-04-25 [2] Bioconductor
# BiocVersion              3.17.1    2022-11-04 [2] Bioconductor
# Biostrings               2.68.1    2023-05-16 [2] Bioconductor
# bit                      4.0.5     2022-11-15 [2] CRAN (R 4.3.1)
# bit64                    4.0.5     2020-08-30 [2] CRAN (R 4.3.1)
# bitops                   1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# blob                     1.2.4     2023-03-17 [2] CRAN (R 4.3.1)
# bslib                    0.5.0     2023-06-09 [2] CRAN (R 4.3.1)
# cachem                   1.0.8     2023-05-01 [2] CRAN (R 4.3.1)
# cli                      3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# codetools                0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
# colorspace               2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# config                   0.3.1     2020-12-17 [2] CRAN (R 4.3.1)
# cowplot                  1.1.1     2020-12-30 [2] CRAN (R 4.3.1)
# crayon                   1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# curl                     5.0.1     2023-06-07 [2] CRAN (R 4.3.1)
# data.table               1.14.8    2023-02-17 [2] CRAN (R 4.3.1)
# DBI                      1.1.3     2022-06-18 [2] CRAN (R 4.3.1)
# dbplyr                   2.3.3     2023-07-07 [2] CRAN (R 4.3.1)
# DelayedArray             0.26.6    2023-07-02 [2] Bioconductor
# DelayedMatrixStats       1.22.1    2023-06-09 [2] Bioconductor
# digest                   0.6.33    2023-07-07 [2] CRAN (R 4.3.1)
# doParallel               1.0.17    2022-02-07 [2] CRAN (R 4.3.1)
# dotCall64                1.0-2     2022-10-03 [2] CRAN (R 4.3.1)
# dplyr                    1.1.2     2023-04-20 [2] CRAN (R 4.3.1)
# dqrng                    0.3.0     2021-05-01 [2] CRAN (R 4.3.1)
# DropletUtils             1.20.0    2023-04-25 [2] Bioconductor
# DT                       0.28      2023-05-18 [2] CRAN (R 4.3.1)
# edgeR                    3.42.4    2023-05-31 [2] Bioconductor
# ellipsis                 0.3.2     2021-04-29 [2] CRAN (R 4.3.1)
# ExperimentHub            2.8.1     2023-07-12 [2] Bioconductor
# fansi                    1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
# fastmap                  1.1.1     2023-02-24 [2] CRAN (R 4.3.1)
# fields                   14.1      2022-08-12 [2] CRAN (R 4.3.1)
# filelock                 1.0.2     2018-10-05 [2] CRAN (R 4.3.1)
# foreach                  1.5.2     2022-02-02 [2] CRAN (R 4.3.1)
# generics                 0.1.3     2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb           * 1.36.1    2023-06-21 [2] Bioconductor
# GenomeInfoDbData         1.2.10    2023-07-20 [2] Bioconductor
# GenomicAlignments        1.36.0    2023-04-25 [2] Bioconductor
# GenomicRanges          * 1.52.0    2023-04-25 [2] Bioconductor
# ggbeeswarm               0.7.2     2023-04-29 [2] CRAN (R 4.3.1)
# ggplot2                  3.4.2     2023-04-03 [2] CRAN (R 4.3.1)
# ggrepel                  0.9.3     2023-02-03 [2] CRAN (R 4.3.1)
# glue                     1.6.2     2022-02-24 [2] CRAN (R 4.3.1)
# golem                    0.4.1     2023-06-05 [2] CRAN (R 4.3.1)
# gridExtra                2.3       2017-09-09 [2] CRAN (R 4.3.1)
# gtable                   0.3.3     2023-03-21 [2] CRAN (R 4.3.1)
# HDF5Array                1.28.1    2023-05-01 [2] Bioconductor
# here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# htmltools                0.5.5     2023-03-23 [2] CRAN (R 4.3.1)
# htmlwidgets              1.6.2     2023-03-17 [2] CRAN (R 4.3.1)
# httpuv                   1.6.11    2023-05-11 [2] CRAN (R 4.3.1)
# httr                     1.4.6     2023-05-08 [2] CRAN (R 4.3.1)
# interactiveDisplayBase   1.38.0    2023-04-25 [2] Bioconductor
# IRanges                * 2.34.1    2023-06-22 [2] Bioconductor
# irlba                    2.3.5.1   2022-10-03 [2] CRAN (R 4.3.1)
# iterators                1.0.14    2022-02-05 [2] CRAN (R 4.3.1)
# jquerylib                0.1.4     2021-04-26 [2] CRAN (R 4.3.1)
# jsonlite                 1.8.7     2023-06-29 [2] CRAN (R 4.3.1)
# KEGGREST                 1.40.0    2023-04-25 [2] Bioconductor
# later                    1.3.1     2023-05-02 [2] CRAN (R 4.3.1)
# lattice                  0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# lazyeval                 0.2.2     2019-03-15 [2] CRAN (R 4.3.1)
# lifecycle                1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
# limma                    3.56.2    2023-06-04 [2] Bioconductor
# lobstr                 * 1.1.2     2022-06-22 [2] CRAN (R 4.3.1)
# locfit                   1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
# magick                   2.7.4     2023-03-09 [2] CRAN (R 4.3.1)
# magrittr                 2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
# maps                     3.4.1     2022-10-30 [2] CRAN (R 4.3.1)
# Matrix                   1.6-0     2023-07-08 [3] CRAN (R 4.3.1)
# MatrixGenerics         * 1.12.2    2023-06-09 [2] Bioconductor
# matrixStats            * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# memoise                  2.0.1     2021-11-26 [2] CRAN (R 4.3.1)
# mime                     0.12      2021-09-28 [2] CRAN (R 4.3.1)
# munsell                  0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
# paletteer                1.5.0     2022-10-19 [2] CRAN (R 4.3.1)
# pillar                   1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# plotly                   4.10.2    2023-06-03 [2] CRAN (R 4.3.1)
# png                      0.1-8     2022-11-29 [2] CRAN (R 4.3.1)
# prettyunits              1.1.1     2020-01-24 [2] CRAN (R 4.3.1)
# promises                 1.2.0.1   2021-02-11 [2] CRAN (R 4.3.1)
# purrr                    1.0.1     2023-01-10 [2] CRAN (R 4.3.1)
# R.methodsS3              1.8.2     2022-06-13 [2] CRAN (R 4.3.1)
# R.oo                     1.25.0    2022-06-12 [2] CRAN (R 4.3.1)
# R.utils                  2.12.2    2022-11-11 [2] CRAN (R 4.3.1)
# R6                       2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# rappdirs                 0.3.3     2021-01-31 [2] CRAN (R 4.3.1)
# RColorBrewer             1.1-3     2022-04-03 [2] CRAN (R 4.3.1)
# Rcpp                     1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RCurl                    1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# rematch2                 2.1.2     2020-05-01 [2] CRAN (R 4.3.1)
# restfulr                 0.0.15    2022-06-16 [2] CRAN (R 4.3.1)
# rhdf5                    2.44.0    2023-04-25 [2] Bioconductor
# rhdf5filters             1.12.1    2023-04-30 [2] Bioconductor
# Rhdf5lib                 1.22.0    2023-04-25 [2] Bioconductor
# rjson                    0.2.21    2022-01-09 [2] CRAN (R 4.3.1)
# rlang                    1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot                2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# Rsamtools                2.16.0    2023-04-25 [2] Bioconductor
# RSQLite                  2.3.1     2023-04-03 [2] CRAN (R 4.3.1)
# rsvd                     1.0.5     2021-04-16 [2] CRAN (R 4.3.1)
# rtracklayer            * 1.60.0    2023-04-25 [2] Bioconductor
# S4Arrays                 1.0.4     2023-05-14 [2] Bioconductor
# S4Vectors              * 0.38.1    2023-05-02 [2] Bioconductor
# sass                     0.4.7     2023-07-15 [2] CRAN (R 4.3.1)
# ScaledMatrix             1.8.1     2023-05-03 [2] Bioconductor
# scales                   1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
# scater                   1.28.0    2023-04-25 [2] Bioconductor
# scuttle                  1.10.1    2023-05-02 [2] Bioconductor
# sessioninfo            * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# shiny                    1.7.4.1   2023-07-06 [2] CRAN (R 4.3.1)
# shinyWidgets             0.7.6     2023-01-08 [2] CRAN (R 4.3.1)
# SingleCellExperiment   * 1.22.0    2023-04-25 [2] Bioconductor
# spam                     2.9-1     2022-08-07 [2] CRAN (R 4.3.1)
# sparseMatrixStats        1.12.2    2023-07-02 [2] Bioconductor
# SpatialExperiment      * 1.10.0    2023-04-25 [2] Bioconductor
# spatialLIBD            * 1.12.0    2023-04-27 [2] Bioconductor
# statmod                  1.5.0     2023-01-06 [2] CRAN (R 4.3.1)
# SummarizedExperiment   * 1.30.2    2023-06-06 [2] Bioconductor
# tibble                   3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
# tidyr                    1.3.0     2023-01-24 [2] CRAN (R 4.3.1)
# tidyselect               1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
# utf8                     1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
# vctrs                    0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
# vipor                    0.4.5     2017-03-22 [2] CRAN (R 4.3.1)
# viridis                  0.6.3     2023-05-03 [2] CRAN (R 4.3.1)
# viridisLite              0.4.2     2023-05-02 [2] CRAN (R 4.3.1)
# XML                      3.99-0.14 2023-03-19 [2] CRAN (R 4.3.1)
# xtable                   1.8-4     2019-04-21 [2] CRAN (R 4.3.1)
# XVector                  0.40.0    2023-04-25 [2] Bioconductor
# yaml                     2.3.7     2023-01-23 [2] CRAN (R 4.3.1)
# zlibbioc                 1.46.0    2023-04-25 [2] Bioconductor
