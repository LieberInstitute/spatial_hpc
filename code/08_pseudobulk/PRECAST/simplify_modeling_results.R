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
# FALSE:9653      FALSE:10318     FALSE:12199     FALSE:9337
# TRUE :2707      TRUE :2042      TRUE :161       TRUE :3023
#     5               6               7               8
# Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:8149      FALSE:7824      FALSE:6368      FALSE:5363
# TRUE :4211      TRUE :4536      TRUE :5992      TRUE :6997
#     10              11              12              13
# Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:8012      FALSE:5968      FALSE:5252      FALSE:3709
# TRUE :4348      TRUE :6392      TRUE :7108      TRUE :8651
#     14
# Mode :logical
# FALSE:8811
# TRUE :3549

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
# [1] 12360   105

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
#1    58.78945  2.044407e-73  2.257157e-73     2.688083 ENSG00000237491 LINC01409
#2   195.02485 4.426314e-132 7.322881e-132     3.947366 ENSG00000228794 LINC01128
#3    40.30790  7.251050e-58  7.507998e-58     1.782555 ENSG00000187634    SAMD11
#4   575.16152 1.108182e-191 3.258892e-191     5.039527 ENSG00000188976     NOC2L
#5    55.86286  3.278979e-71  3.578332e-71     2.586253 ENSG00000187961    KLHL17
#6  1341.23768 3.925086e-240 2.008032e-239     5.652542 ENSG00000188290      HES4

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
