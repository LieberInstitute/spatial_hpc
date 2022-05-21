setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("scran")
})

load(file = here::here("processed-data", "pilot_data_checks", "spe_raw.Rdata"))

## grid
vis_grid_clus(
    spe = spe_raw,
    clustervar = "in_tissue",
    pdf = here::here("plots", "pilot_data_checks", "in_tissue_grid.pdf"),
    sort_clust = FALSE,
    colors = c("TRUE" = "grey90", "FALSE" = "orange")
)



head(table(spe_raw$sum_umi[which(spatialData(spe_raw)$in_tissue == "FALSE")]))
# 4  8  9 10 12 13
# 1  1  1  2  1  2

## out of tissue
vis_grid_gene(
    spe = spe_raw[, which(spatialData(spe_raw)$in_tissue == "FALSE")],
    geneid = "sum_umi",
    pdf = here::here("plots", "pilot_data_checks", "out_tissue_sum_umi.pdf"),
    assayname = "counts"
)

summary(spe_raw$sum_umi[which(spatialData(spe_raw)$in_tissue == "FALSE")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 4.0   129.0   229.0   385.8   400.0 20096.0

vis_grid_gene(
    spe = spe_raw[, which(spatialData(spe_raw)$in_tissue == "FALSE")],
    geneid = "sum_gene",
    pdf = here::here("plots", "pilot_data_checks", "out_tissue_sum_gene.pdf"),
    assayname = "counts"
)

summary(spe_raw$sum_gene[which(spatialData(spe_raw)$in_tissue == "FALSE")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 4.0   107.0   178.0   265.9   297.0  6222.0

vis_grid_gene(
    spe = spe_raw[, which(spatialData(spe_raw)$in_tissue == "FALSE")],
    geneid = "expr_chrM_ratio",
    pdf = here::here("plots", "pilot_data_checks", "out_tissue_expr_chrM_ratio.pdf"),
    assayname = "counts"
)

summary(spe_raw$expr_chrM_ratio[which(spatialData(spe_raw)$in_tissue == "FALSE")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.1205  0.1531  0.1585  0.1883  0.4800

## in tissue
vis_grid_gene(
    spe = spe_raw[, which(spatialData(spe_raw)$in_tissue == "TRUE")],
    geneid = "sum_umi",
    pdf = here::here("plots", "pilot_data_checks", "in_tissue_sum_umi.pdf"),
    assayname = "counts"
)


vis_grid_gene(
    spe = spe_raw[, which(spatialData(spe_raw)$in_tissue == "TRUE")],
    geneid = "sum_gene",
    pdf = here::here("plots", "pilot_data_checks", "in_tissue_sum_gene.pdf"),
    assayname = "counts"
)


vis_grid_gene(
    spe = spe_raw[, which(spatialData(spe_raw)$in_tissue == "TRUE")],
    geneid = "expr_chrM_ratio",
    pdf = here::here("plots", "pilot_data_checks", "in_tissue_expr_chrM_ratio.pdf"),
    assayname = "counts"
)


## all
vis_grid_gene(
    spe = spe_raw,
    geneid = "sum_umi",
    pdf = here::here("plots", "pilot_data_checks", "all_sum_umi.pdf"),
    assayname = "counts"
)
summary(spe_raw$sum_umi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 4     552    1893    2735    3712   63545


vis_grid_gene(
    spe = spe_raw,
    geneid = "sum_gene",
    pdf = here::here("plots", "pilot_data_checks", "all_sum_gene.pdf"),
    assayname = "counts"
)
summary(spe_raw$sum_gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 4     387    1117    1375    2009    9770

vis_grid_gene(
    spe = spe_raw,
    geneid = "expr_chrM_ratio",
    pdf = here::here("plots", "pilot_data_checks", "all_expr_chrM_ratio.pdf"),
    assayname = "counts"
)

summary(spe_raw$expr_chrM_ratio)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00000 0.09051 0.13168 0.13877 0.17647 0.53391

## Quality control (scran)
qcstats <- perCellQCMetrics(spe_raw, subsets = list(
    Mito = which(seqnames(spe_raw) == "chrM")
))
qcfilter <- quickPerCellQC(qcstats, sub.fields = "subsets_Mito_percent")
colSums(as.matrix(qcfilter))
# low_lib_size            low_n_features high_subsets_Mito_percent
# 598                       782                       489
# discard
# 1233

## Prior to dropping spots with 0 counts and checking for high chrM, this was the output:

spe_raw$scran_discard <-
    factor(qcfilter$discard, levels = c("TRUE", "FALSE"))
spe_raw$scran_low_lib_size <-
    factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
spe_raw$scran_low_n_features <-
    factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
spe_raw$scran_high_subsets_Mito_percent <-
    factor(qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))
save(spe_raw, file = here::here("processed-data", "pilot_data_checks", "spe_final.Rdata"))

for (i in colnames(qcfilter)) {
    vis_grid_clus(
        spe = spe_raw,
        clustervar = paste0("scran_", i),
        pdf = here::here("plots", "pilot_data_checks", paste0("scran_", i, ".pdf")),
        sort_clust = FALSE,
        colors = c("FALSE" = "grey90", "TRUE" = "orange")
    )
}
