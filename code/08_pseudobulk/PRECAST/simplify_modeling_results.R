###############################
# spatial_HPC project
# Parse DE analysis results
# Anthony Ramnauth, Dec 21 2022
###############################

setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(scater)
    library(sessioninfo)
})

# Load modeling results
load(file = here::here("processed-data", "08_pseudobulk", "PRECAST",
    "DE_list_captureArea_adjBrnum.Rdata"))

## load spe data
load(file = here::here("processed-data", "08_pseudobulk", "PRECAST",
    "spe_pseudo_captureArea_wo_9-15-NA_Fncells50.Rdata"))

## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
    x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts) <- rownames(eb_contrasts)
fdrs0_contrasts <- apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the t-stats
t0_contrasts <- sapply(eb0_list, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts) <- rownames(eb_contrasts)
summary(fdrs0_contrasts < 0.05)
#     1               2               3               4
# Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:10083     FALSE:10685     FALSE:12300     FALSE:9517
# TRUE :2271      TRUE :1669      TRUE :54        TRUE :2837
#     5               6               7               8
# Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:8425      FALSE:8175      FALSE:6730      FALSE:5374
# TRUE :3929      TRUE :4179      TRUE :5624      TRUE :6980
#     10              11              12              13
# Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:8536      FALSE:6552      FALSE:5566      FALSE:2233
# TRUE :3818      TRUE :5802      TRUE :6788      TRUE :10121
#     14
# Mode :logical
# FALSE:9234
# TRUE :3120

# Merge statistics
f_merge <- function(p, fdr, t) {
    colnames(p) <- paste0("p_value_", colnames(p))
    colnames(fdr) <- paste0("fdr_", colnames(fdr))
    colnames(t) <- paste0("t_stat_", colnames(t))
    res <- as.data.frame(cbind(t, p, fdr))
    res$ensembl <- rownames(res)
    ## Check it's all in order
    stopifnot(identical(rownames(res), rownames(spe_pseudo)))
    res$gene <- rowData(spe_pseudo)$gene_name
    rownames(res) <- NULL
    return(res)
}

results_specificity <-
    f_merge(p = pvals0_contrasts, fdr = fdrs0_contrasts, t = t0_contrasts)
options(width = 400)
head(results_specificity)

pvals_contrasts <- eb_contrasts$p.value
fdrs_contrasts <- apply(pvals_contrasts, 2, p.adjust, "fdr")
dim(pvals_contrasts)
# [1] 12354   105

summary(fdrs_contrasts < 0.05)

results_pairwise <-
    f_merge(p = pvals_contrasts, fdr = fdrs_contrasts, t = eb_contrasts$t)
colnames(results_pairwise)
sort(colSums(fdrs_contrasts < 0.05))

f_sig <- function(type, cut = 0.05) {
    cbind(
        "n" = addmargins(table(f_stats[[type]] < cut)),
        "ratio" = addmargins(table(f_stats[[type]] < cut)) / nrow(f_stats)
    )
}
f_sig("full_fdr")

# Match the colnames to the new style
f_rename <- function(x, old, new = old) {
    old_patt <- paste0("_", old, "$")
    i <- grep(old_patt, colnames(x))
    tmp <- gsub(old_patt, "", colnames(x)[i])
    tmp <- paste0(new, "_", tmp)
    colnames(x)[i] <- tmp
    return(x)
}

results_anova <-
    f_rename(f_rename(f_rename(
        f_rename(f_stats, "f", "f_stat"), "p_value"
    ), "fdr"), "Amean")
head(results_anova)
#  f_stat_full  p_value_full      fdr_full full_AveExpr         ensembl      gene
#1    60.95011  5.602754e-76  6.196635e-76     2.636586 ENSG00000237491 LINC01409
#2   174.57190 2.278567e-128 3.645825e-128     3.879536 ENSG00000228794 LINC01128
#3    43.82572  7.093350e-62  7.383825e-62     1.762849 ENSG00000187634    SAMD11
#4   535.97863 2.236853e-191 6.717083e-191     5.026624 ENSG00000188976     NOC2L
#5    56.53032  1.217997e-72  1.322941e-72     2.537504 ENSG00000187961    KLHL17
#6  1523.38833 9.685692e-253 5.944214e-252     5.640188 ENSG00000188290      HES4

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_specificity,
    "pairwise" = results_pairwise
)

# Save modeling results
saveRDS(modeling_results, file = here::here("processed-data", "08_pseudobulk", "PRECAST",
    "modeling_results.rds"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
