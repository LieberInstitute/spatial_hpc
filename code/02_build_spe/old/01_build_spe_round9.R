setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc")
# cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/
suppressPackageStartupMessages(library("here"))
# remotes::install_github("drighelli/SpatialExperiment")
# remotes::install_github("LieberInstitute/spatialLIBD")
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("lobstr"))
suppressPackageStartupMessages(library("sessioninfo"))

## Create output directories
dir_rdata <- here::here("processed-data", "02_build_spe")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Define some info for the samples
sample_info <- data.frame(
    sample_id = c(
        "Br2720_A1",
        "Br2720_B1",
        "Br2720_C1",
        "Br2720_D1"
    )
)
sample_info$subject <- sample_info$sample_id
sample_info$sample_path <-
    file.path(
        here::here("processed-data", "01_spaceranger", "spaceranger_2023-01-31_round9"),
        sample_info$sample_id,
        "outs"
    )
stopifnot(all(file.exists(sample_info$sample_path)))

## Define the donor info using information from
## /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/raw-data/sample_info/demographicInfo_Geo.csv
donor_info <- data.frame(
    subject = c("Br2720_A1","Br2720_B1", "Br2720_C1", "Br2720_D1"),
    age = c(48.2, 48.2, 48.2, 48.2),
    sex = c("F", "F", "F", "F"),
    race = c("EA/CAUC", "EA/CAUC", "EA/CAUC", "EA/CAUC"),
    diagnosis = c("Control", "Control", "Control", "Control"),
    rin = c(7.5, 7.5, 7.5, 7.5),
    pmi = c(25.5, 25.5, 25.5, 25.5)
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

# 2023-02-13 14:28:52 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
# 2023-02-13 14:29:23 read10xVisiumAnalysis: reading analysis output from SpaceRanger
# 2023-02-13 14:29:25 add10xVisiumAnalysis: adding analysis output from SpaceRanger
# 2023-02-13 14:29:26 rtracklayer::import: reading the reference GTF file
# 2023-02-13 14:30:03 adding gene information to the SPE object
# 2023-02-13 14:30:03 adding information used by spatialLIBD
# "2023-02-13 14:31:25 EST"

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
                "01_spaceranger",
                "spaceranger_2023-01-31_round9",
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

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)
length(no_expr)
# [1] 11473
length(no_expr) / nrow(spe) * 100
# [1] 31.34614
spe <- spe[-no_expr, ]


## For visualizing this later with spatialLIBD
spe$overlaps_tissue <-
    factor(ifelse(spe$in_tissue, "in", "out"))

## Save with and without dropping spots outside of the tissue
spe_raw_round9 <- spe

saveRDS(spe_raw, file.path(dir_rdata, "spe_raw_round9.rds"))

## Size in Gb
lobstr::obj_size(spe_raw_round9)
# 668.28 MB


## Now drop the spots outside the tissue
spe <- spe_raw_round9[, spe_raw_round9$in_tissue]
dim(spe)
# [1] 25128 15641
## Remove spots without counts
if (any(colSums(counts(spe)) == 0)) {
    message("removing spots without counts for spe")
    spe <- spe[, -which(colSums(counts(spe)) == 0)]
    dim(spe)
}

lobstr::obj_size(spe)
# 659.29 MB

saveRDS(spe, file.path(dir_rdata, "spe_round9.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2023-02-13 14:44:20 EST"
#     user   system  elapsed
#  283.754    5.259 1547.371
# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.2.1 Patched (2022-08-30 r82775)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2023-02-13
#  pandoc   2.19.2 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2/bin/pandoc
#
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package                * version   date (UTC) lib source
#  AnnotationDbi            1.58.0    2022-04-26 [2] Bioconductor
#  AnnotationHub            3.4.0     2022-04-26 [2] Bioconductor
#  assertthat               0.2.1     2019-03-21 [2] CRAN (R 4.2.1)
#  attempt                  0.3.1     2020-05-03 [2] CRAN (R 4.2.1)
#  beachmat                 2.12.0    2022-04-26 [2] Bioconductor
#  beeswarm                 0.4.0     2021-06-01 [2] CRAN (R 4.2.1)
#  benchmarkme              1.0.8     2022-06-12 [2] CRAN (R 4.2.1)
#  benchmarkmeData          1.0.4     2020-04-23 [2] CRAN (R 4.2.1)
#  Biobase                * 2.56.0    2022-04-26 [2] Bioconductor
#  BiocFileCache            2.4.0     2022-04-26 [2] Bioconductor
#  BiocGenerics           * 0.42.0    2022-04-26 [2] Bioconductor
#  BiocIO                   1.6.0     2022-04-26 [2] Bioconductor
#  BiocManager              1.30.18   2022-05-18 [2] CRAN (R 4.2.1)
#  BiocNeighbors            1.14.0    2022-04-26 [2] Bioconductor
#  BiocParallel             1.30.3    2022-06-05 [2] Bioconductor
#  BiocSingular             1.12.0    2022-04-26 [2] Bioconductor
#  BiocVersion              3.15.2    2022-03-29 [2] Bioconductor
#  Biostrings               2.64.1    2022-08-18 [2] Bioconductor
#  bit                      4.0.4     2020-08-04 [2] CRAN (R 4.2.1)
#  bit64                    4.0.5     2020-08-30 [2] CRAN (R 4.2.1)
#  bitops                   1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
#  blob                     1.2.3     2022-04-10 [2] CRAN (R 4.2.1)
#  bslib                    0.4.0     2022-07-16 [2] CRAN (R 4.2.1)
#  cachem                   1.0.6     2021-08-19 [2] CRAN (R 4.2.1)
#  cli                      3.3.0     2022-04-25 [2] CRAN (R 4.2.1)
#  codetools                0.2-18    2020-11-04 [3] CRAN (R 4.2.1)
#  colorspace               2.0-3     2022-02-21 [2] CRAN (R 4.2.1)
#  config                   0.3.1     2020-12-17 [2] CRAN (R 4.2.1)
#  cowplot                  1.1.1     2020-12-30 [2] CRAN (R 4.2.1)
#  crayon                   1.5.1     2022-03-26 [2] CRAN (R 4.2.1)
#  curl                     4.3.2     2021-06-23 [2] CRAN (R 4.2.1)
#  data.table               1.14.2    2021-09-27 [2] CRAN (R 4.2.1)
#  DBI                      1.1.3     2022-06-18 [2] CRAN (R 4.2.1)
#  dbplyr                   2.2.1     2022-06-27 [2] CRAN (R 4.2.1)
#  DelayedArray             0.22.0    2022-04-26 [2] Bioconductor
#  DelayedMatrixStats       1.18.0    2022-04-26 [2] Bioconductor
#  desc                     1.4.1     2022-03-06 [2] CRAN (R 4.2.1)
#  digest                   0.6.29    2021-12-01 [2] CRAN (R 4.2.1)
#  doParallel               1.0.17    2022-02-07 [2] CRAN (R 4.2.1)
#  dotCall64                1.0-1     2021-02-11 [2] CRAN (R 4.2.1)
#  dplyr                    1.0.9     2022-04-28 [2] CRAN (R 4.2.1)
#  dqrng                    0.3.0     2021-05-01 [2] CRAN (R 4.2.1)
#  DropletUtils             1.16.0    2022-04-26 [2] Bioconductor
#  DT                       0.24      2022-08-09 [2] CRAN (R 4.2.1)
#  edgeR                    3.38.4    2022-08-07 [2] Bioconductor
#  ellipsis                 0.3.2     2021-04-29 [2] CRAN (R 4.2.1)
#  ExperimentHub            2.4.0     2022-04-26 [2] Bioconductor
#  fansi                    1.0.3     2022-03-24 [2] CRAN (R 4.2.1)
#  fastmap                  1.1.0     2021-01-25 [2] CRAN (R 4.2.1)
#  fields                   14.1      2022-08-12 [2] CRAN (R 4.2.1)
#  filelock                 1.0.2     2018-10-05 [2] CRAN (R 4.2.1)
#  foreach                  1.5.2     2022-02-02 [2] CRAN (R 4.2.1)
#  fs                       1.5.2     2021-12-08 [2] CRAN (R 4.2.1)
#  generics                 0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
#  GenomeInfoDb           * 1.32.3    2022-08-09 [2] Bioconductor
#  GenomeInfoDbData         1.2.8     2022-08-30 [2] Bioconductor
#  GenomicAlignments        1.32.1    2022-07-24 [2] Bioconductor
#  GenomicRanges          * 1.48.0    2022-04-26 [2] Bioconductor
#  ggbeeswarm               0.6.0     2017-08-07 [2] CRAN (R 4.2.1)
#  ggplot2                  3.3.6     2022-05-03 [2] CRAN (R 4.2.1)
#  ggrepel                  0.9.1     2021-01-15 [2] CRAN (R 4.2.1)
#  glue                     1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
#  golem                    0.3.3     2022-07-13 [2] CRAN (R 4.2.1)
#  gridExtra                2.3       2017-09-09 [2] CRAN (R 4.2.1)
#  gtable                   0.3.0     2019-03-25 [2] CRAN (R 4.2.1)
#  HDF5Array                1.24.2    2022-08-02 [2] Bioconductor
#  here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
#  htmltools                0.5.3     2022-07-18 [2] CRAN (R 4.2.1)
#  htmlwidgets              1.5.4     2021-09-08 [2] CRAN (R 4.2.1)
#  httpuv                   1.6.5     2022-01-05 [2] CRAN (R 4.2.1)
#  httr                     1.4.4     2022-08-17 [2] CRAN (R 4.2.1)
#  interactiveDisplayBase   1.34.0    2022-04-26 [2] Bioconductor
#  IRanges                * 2.30.1    2022-08-18 [2] Bioconductor
#  irlba                    2.3.5     2021-12-06 [2] CRAN (R 4.2.1)
#  iterators                1.0.14    2022-02-05 [2] CRAN (R 4.2.1)
#  jquerylib                0.1.4     2021-04-26 [2] CRAN (R 4.2.1)
#  jsonlite                 1.8.0     2022-02-22 [2] CRAN (R 4.2.1)
#  KEGGREST                 1.36.3    2022-07-12 [2] Bioconductor
#  knitr                    1.40      2022-08-24 [2] CRAN (R 4.2.1)
#  later                    1.3.0     2021-08-18 [2] CRAN (R 4.2.1)
#  lattice                  0.20-45   2021-09-22 [3] CRAN (R 4.2.1)
#  lazyeval                 0.2.2     2019-03-15 [2] CRAN (R 4.2.1)
#  lifecycle                1.0.1     2021-09-24 [2] CRAN (R 4.2.1)
#  limma                    3.52.2    2022-06-19 [2] Bioconductor
#  lobstr                 * 1.1.2     2022-06-22 [2] CRAN (R 4.2.1)
#  locfit                   1.5-9.6   2022-07-11 [2] CRAN (R 4.2.1)
#  magick                   2.7.3     2021-08-18 [2] CRAN (R 4.2.1)
#  magrittr                 2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
#  maps                     3.4.0     2021-09-25 [2] CRAN (R 4.2.1)
#  Matrix                   1.4-1     2022-03-23 [3] CRAN (R 4.2.1)
#  MatrixGenerics         * 1.8.1     2022-06-26 [2] Bioconductor
#  matrixStats            * 0.62.0    2022-04-19 [2] CRAN (R 4.2.1)
#  memoise                  2.0.1     2021-11-26 [2] CRAN (R 4.2.1)
#  mime                     0.12      2021-09-28 [2] CRAN (R 4.2.1)
#  munsell                  0.5.0     2018-06-12 [2] CRAN (R 4.2.1)
#  paletteer                1.4.1     2022-08-15 [2] CRAN (R 4.2.1)
#  pillar                   1.8.1     2022-08-19 [2] CRAN (R 4.2.1)
#  pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
#  pkgload                  1.3.0     2022-06-27 [2] CRAN (R 4.2.1)
#  plotly                   4.10.0    2021-10-09 [2] CRAN (R 4.2.1)
#  png                      0.1-7     2013-12-03 [2] CRAN (R 4.2.1)
#  prettyunits              1.1.1     2020-01-24 [2] CRAN (R 4.2.1)
#  promises                 1.2.0.1   2021-02-11 [2] CRAN (R 4.2.1)
#  purrr                    0.3.4     2020-04-17 [2] CRAN (R 4.2.1)
#  R.methodsS3              1.8.2     2022-06-13 [2] CRAN (R 4.2.1)
#  R.oo                     1.25.0    2022-06-12 [2] CRAN (R 4.2.1)
#  R.utils                  2.12.0    2022-06-28 [2] CRAN (R 4.2.1)
#  R6                       2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
#  rappdirs                 0.3.3     2021-01-31 [2] CRAN (R 4.2.1)
#  RColorBrewer             1.1-3     2022-04-03 [2] CRAN (R 4.2.1)
#  Rcpp                     1.0.9     2022-07-08 [2] CRAN (R 4.2.1)
#  RCurl                    1.98-1.8  2022-07-30 [2] CRAN (R 4.2.1)
#  rematch2                 2.1.2     2020-05-01 [2] CRAN (R 4.2.1)
#  restfulr                 0.0.15    2022-06-16 [2] CRAN (R 4.2.1)
#  rhdf5                    2.40.0    2022-04-26 [2] Bioconductor
#  rhdf5filters             1.8.0     2022-04-26 [2] Bioconductor
#  Rhdf5lib                 1.18.2    2022-05-15 [2] Bioconductor
#  rjson                    0.2.21    2022-01-09 [2] CRAN (R 4.2.1)
#  rlang                    1.0.4     2022-07-12 [2] CRAN (R 4.2.1)
#  roxygen2                 7.2.1     2022-07-18 [2] CRAN (R 4.2.1)
#  rprojroot                2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
#  Rsamtools                2.12.0    2022-04-26 [2] Bioconductor
#  RSQLite                  2.2.16    2022-08-17 [2] CRAN (R 4.2.1)
#  rstudioapi               0.14      2022-08-22 [2] CRAN (R 4.2.1)
#  rsvd                     1.0.5     2021-04-16 [2] CRAN (R 4.2.1)
#  rtracklayer            * 1.56.1    2022-06-23 [2] Bioconductor
#  S4Vectors              * 0.34.0    2022-04-26 [2] Bioconductor
#  sass                     0.4.2     2022-07-16 [2] CRAN (R 4.2.1)
#  ScaledMatrix             1.4.0     2022-04-26 [2] Bioconductor
#  scales                   1.2.1     2022-08-20 [2] CRAN (R 4.2.1)
#  scater                   1.24.0    2022-04-26 [2] Bioconductor
#  scuttle                  1.6.3     2022-08-23 [2] Bioconductor
#  sessioninfo            * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
#  shiny                    1.7.2     2022-07-19 [2] CRAN (R 4.2.1)
#  shinyWidgets             0.7.2     2022-08-07 [2] CRAN (R 4.2.1)
#  SingleCellExperiment   * 1.18.0    2022-04-26 [2] Bioconductor
#  spam                     2.9-1     2022-08-07 [2] CRAN (R 4.2.1)
#  sparseMatrixStats        1.8.0     2022-04-26 [2] Bioconductor
#  SpatialExperiment      * 1.6.1     2022-08-09 [2] Bioconductor
#  spatialLIBD            * 1.8.10    2022-07-26 [2] Bioconductor
#  stringi                  1.7.8     2022-07-11 [2] CRAN (R 4.2.1)
#  stringr                  1.4.1     2022-08-20 [2] CRAN (R 4.2.1)
#  SummarizedExperiment   * 1.26.1    2022-04-29 [2] Bioconductor
#  tibble                   3.1.8     2022-07-22 [2] CRAN (R 4.2.1)
#  tidyr                    1.2.0     2022-02-01 [2] CRAN (R 4.2.1)
#  tidyselect               1.1.2     2022-02-21 [2] CRAN (R 4.2.1)
#  usethis                  2.1.6     2022-05-25 [2] CRAN (R 4.2.1)
#  utf8                     1.2.2     2021-07-24 [2] CRAN (R 4.2.1)
#  vctrs                    0.4.1     2022-04-13 [2] CRAN (R 4.2.1)
#  vipor                    0.4.5     2017-03-22 [2] CRAN (R 4.2.1)
#  viridis                  0.6.2     2021-10-13 [2] CRAN (R 4.2.1)
#  viridisLite              0.4.1     2022-08-22 [2] CRAN (R 4.2.1)
#  xfun                     0.32      2022-08-10 [2] CRAN (R 4.2.1)
#  XML                      3.99-0.10 2022-06-09 [2] CRAN (R 4.2.1)
#  xml2                     1.3.3     2021-11-30 [2] CRAN (R 4.2.1)
#  xtable                   1.8-4     2019-04-21 [2] CRAN (R 4.2.1)
#  XVector                  0.36.0    2022-04-26 [2] Bioconductor
#  yaml                     2.3.5     2022-02-21 [2] CRAN (R 4.2.1)
#  zlibbioc                 1.42.0    2022-04-26 [2] Bioconductor
#
#  [1] /users/hdivecha/R/4.2
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2/R/4.2/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2/R/4.2/lib64/R/library

