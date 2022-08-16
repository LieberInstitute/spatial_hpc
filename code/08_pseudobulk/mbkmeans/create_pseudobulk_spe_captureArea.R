setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(spatialLIBD)
  library(here)
  library(edgeR)
  library(scuttle)
  library(scater)
  library(scran)
  library(dplyr)
  library(PCAtools)
  library(sessioninfo)
  library(gridExtra)
})

# Create directory for BayesSpace pseudo-bulked spe object
dir_rdata <- here::here("processed-data", "08_pseudobulk", "mbkmeans")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

# Create directory for pseudobulked plots
dir_plots <- here::here("plots", "08_pseudobulk", "mbkmeans")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
load(file = here::here("processed-data", "06_Clustering", "spe_modify.Rdata"))
load(file = here("processed-data", "06_Clustering", "mbkmeans.Rdata"))

spe$kmeans <- km_res[[13]]$Clusters
## Pseudo-bulk for BayesSpace k = 11 results
spe_pseudo <- aggregateAcrossCells(
  spe,
  DataFrame(
    mbkmeans = colData(spe)$kmeans,
    sample_id = spe$sample_id
  ))

spe_pseudo$mbkmeans <- factor(spe_pseudo$mbkmeans)

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr <- filterByExpr(spe_pseudo)
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$mbkmeans)

summary(rowData(spe_pseudo)$high_expr)
summary(rowData(spe_pseudo)$high_expr_group_sample_id)
summary(rowData(spe_pseudo)$high_expr_group_cluster)

with(rowData(spe_pseudo), table(high_expr_group_sample_id, high_expr_group_cluster))

## Now filter
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_cluster, ]
dim(spe_pseudo)

# Store the log normalized counts on the spe object
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)

# Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo)))

# Fix the column names. DGEList will have samples name as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)

# Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo) <- x
dim(spe_pseudo)

rm(x)

#run PCA
pca <- prcomp(t(assays(spe_pseudo)$logcounts))

message(Sys.time(), " % of variance explained for the top 20 PCs:")
metadata(spe_pseudo)
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
metadata(spe_pseudo)
# $PCA_var_explained
# [1] 43.400  8.560  2.630  2.010  1.570  1.220  0.962  0.786  0.639  0.602
# [11]  0.487  0.448  0.443  0.438  0.429  0.425  0.416  0.395  0.389  0.384

pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 43.400  8.560  2.630  2.010  1.570  1.220  0.962  0.786  0.639  0.602
# [11]  0.487  0.448  0.443  0.438  0.429  0.425  0.416  0.395  0.389  0.384
# [21]  0.376  0.371  0.367  0.355  0.347  0.342  0.340  0.329  0.327  0.323
# [31]  0.316  0.309  0.307  0.302  0.299  0.290  0.286  0.282  0.276  0.273
# [41]  0.260  0.256  0.251  0.248  0.243  0.239  0.236  0.234  0.232  0.226

# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "pseudobulk_captureArea_PCA_mbkmenas17.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "brnum", ncomponents = 12, point_size = 3, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "mbkmeans", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
dev.off()

pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "pseudobulk_captureArea_PCA_mbkmenas17_2.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "brnum", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "mbkmeans", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
dev.off()

pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "pseudobulk_captureArea_PCA_mbkmeans17_4.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "brnum", ncomponents = 4, point_size = 4, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "mbkmeans", ncomponents = 4, point_size = 4, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 4, point_size = 4, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 4, point_size = 4, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 4, point_size = 4, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)

dev.off()

## Compute some reduced dims
set.seed(20220423)
spe_pseudo <- scater::runMDS(spe_pseudo, ncomponents = 20)
spe_pseudo <- scater::runPCA(spe_pseudo, name = "runPCA")

####plot explanatory variables ####

#uses linear regression model
vars <- getVarianceExplained(spe_pseudo, variables=c("brnum", "mbkmeans","sample_id","age","sex"))
head(vars)
# brnum mbkmeans sample_id         age         sex
# ENSG00000237491 5.1131169 26.86047  9.033304 0.076173031 0.049237682
# ENSG00000228794 0.9484676 45.58115  3.378215 0.007751482 0.153853642
# ENSG00000187634 5.1216385 32.67619 12.519473 0.299421815 2.886575179
# ENSG00000188976 1.6498779 60.27809  5.290631 0.071596550 0.002581295
# ENSG00000187961 1.4295300 26.71936  5.246689 0.037105746 0.003321983
# ENSG00000188290 6.2041552 60.43831 10.528557 0.038585217 4.018723476

pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "plot_explanatory_vars_captureArea_mbkmeans17.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save file
save(spe_pseudo, file = here::here("processed-data", "08_pseudobulk", "mbkmeans", "spe_pseudo_captureArea_mbkmeans17.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()