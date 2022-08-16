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
table(spe$sample_id,spe$kmeans)

speb = spe[, which(spe$kmeans != "4")]
speb = speb[, which(speb$kmeans != "12")]
table(speb$sample_id,speb$kmeans)
## Pseudo-bulk for BayesSpace k = 11 results
spe_pseudo <- aggregateAcrossCells(
  speb,
  DataFrame(
    mbkmeans = colData(speb)$kmeans,
    sample_id = speb$sample_id
  ))

spe_pseudo$mbkmeans <- factor(spe_pseudo$mbkmeans)
table(spe_pseudo$sample_id,spe_pseudo$mbkmeans)
table(spe_pseudo$sample_id,spe_pseudo$ncells)
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 10]

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$mbkmeans)

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
# [1] 17.500 10.500  4.290  2.740  2.680  1.920  1.440  1.220  1.120  0.965
# [11]  0.958  0.914  0.889  0.859  0.828  0.790  0.764  0.746  0.737  0.732

pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 17.500 10.500  4.290  2.740  2.680  1.920  1.440  1.220  1.120  0.965
# [11]  0.958  0.914  0.889  0.859  0.828  0.790  0.764  0.746  0.737  0.732
# [21]  0.710  0.694  0.691  0.675  0.655  0.644  0.633  0.613  0.592  0.577
# [31]  0.571  0.569  0.555  0.546  0.518  0.506  0.491  0.473  0.471  0.448
# [41]  0.439  0.439  0.431  0.425  0.417  0.408  0.389  0.376  0.369  0.362

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
# brnum mbkmeans sample_id          age         sex
# ENSG00000228794  1.4221316 22.53474  3.999842 2.698858e-01 0.327654891
# ENSG00000188976  0.9836954 34.10381  5.591324 1.017375e-01 0.001846929
# ENSG00000188290 13.7472569 31.99271 21.244127 3.792593e-01 6.566212901
# ENSG00000187608 17.1798913 22.14732 25.604535 1.113442e-05 0.677990670
# ENSG00000188157 10.8343399 17.57840 22.753861 2.980362e-01 0.019242937
# ENSG00000131591  3.6516595 13.35010 11.605515 1.634773e-01 1.139441196

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