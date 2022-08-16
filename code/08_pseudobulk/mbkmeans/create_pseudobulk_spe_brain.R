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
  library(SingleCellExperiment)
})


# Load SPE
load(file = here::here("processed-data", "06_Clustering", "spe_modify.Rdata"))
load(file = here("processed-data", "06_Clustering", "mbkmeans.Rdata"))

spe$kmeans <- km_res[[13]]$Clusters
## Pseudo-bulk for mbkmeans k = 17 results
sce <- as(spe, "SingleCellExperiment")
spe_pseudo <- aggregateAcrossCells(
  sce,
  id=DataFrame(
    mbkmeans = colData(sce)$kmeans,
    brnum = colData(sce)$brnum)
)

spe_pseudo$mbkmeans <- factor(spe_pseudo$mbkmeans)
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 11]
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
#[1] 68.100  5.680  2.010  1.490  1.290  1.170  1.050  0.978  0.936  0.896
#[11]  0.853  0.799  0.725  0.674  0.653  0.552  0.500  0.466  0.440  0.408

pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 68.100  5.680  2.010  1.490  1.290  1.170  1.050  0.978  0.936  0.896
# [11]  0.853  0.799  0.725  0.674  0.653  0.552  0.500  0.466  0.440  0.408
# [21]  0.373  0.338  0.324  0.313  0.309  0.308  0.287  0.285  0.280  0.267
# [31]  0.260  0.248  0.241  0.226  0.220  0.216  0.204  0.196  0.190  0.187
# [41]  0.174  0.162  0.156  0.149  0.137  0.119  0.115  0.112  0.105  0.104

# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "pseudobulk_brain_PCA_mbkmenas17.pdf"), width = 14, height = 14)
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

pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "pseudobulk_brain_PCA_mbkmenas17_2.pdf"), width = 14, height = 14)
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

pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "pseudobulk_brain_PCA_mbkmeans17_4.pdf"), width = 14, height = 14)
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

# brnum mbkmeans sample_id         age       sex
# ENSG00000241860 21.073382 47.30329       NaN 0.005758258 2.5406981
# ENSG00000237491  6.826028 60.52541       NaN 0.001834419 0.8945204
# ENSG00000228794  1.464072 93.27664       NaN 0.101931629 0.1325683
# ENSG00000230368 21.690064 33.00999       NaN 2.809853149 0.6668615
# ENSG00000223764  4.135560 52.18557       NaN 0.243995288 1.4884063
# ENSG00000187634  7.494925 65.15615       NaN 0.639563646 2.8006216

pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "plot_explanatory_vars_brain_mbkmeans17.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save file
save(spe_pseudo, file = here::here("processed-data", "08_pseudobulk", "mbkmeans", "spe_pseudo_brain_mbkmeans17.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()