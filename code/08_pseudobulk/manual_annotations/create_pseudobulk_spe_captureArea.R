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
dir_rdata <- here::here("processed-data", "08_pseudobulk", "manual_annotations")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

# Create directory for pseudobulked plots
dir_plots <- here::here("plots", "08_pseudobulk", "manual_annotations")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
load(file = here::here("processed-data", "06_Clustering", "spe_modify.Rdata"))

table(spe$sample_id,spe$ManualAnnotations)
## Pseudo-bulk for manualannotations k = 15 results
spe_pseudo <- aggregateAcrossCells(
  spe,
  DataFrame(
    manual_annotations = colData(spe)$ManualAnnotation,
    sample_id = spe$sample_id
  ))

spe_pseudo$manual_annotations <- factor(spe_pseudo$manual_annotations)
table(spe_pseudo$sample_id,spe_pseudo$manual_annotations)
table(spe_pseudo$sample_id,spe_pseudo$ncells)
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 10]

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$manual_annotations)

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
# [1] 12.300 11.600  4.530  3.250  2.590  2.290  1.780  1.660  1.650  1.530
# [11]  1.440  1.390  1.320  1.220  1.130  1.080  1.050  0.987  0.968  0.914

pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 12.300 11.600  4.530  3.250  2.590  2.290  1.780  1.660  1.650  1.530
# [11]  1.440  1.390  1.320  1.220  1.130  1.080  1.050  0.987  0.968  0.914
# [21]  0.845  0.805  0.765  0.753  0.709  0.683  0.657  0.638  0.633  0.627
# [31]  0.618  0.588  0.571  0.559  0.554  0.542  0.518  0.509  0.503  0.492
# [41]  0.479  0.475  0.465  0.455  0.453  0.448  0.433  0.430  0.429  0.419

# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "pseudobulk_captureArea_PCA_manual_annotations.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "brnum", ncomponents = 12, point_size = 3, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "manual_annotations", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
dev.off()

pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "pseudobulk_captureArea_PCA_manual_annotations_2.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "brnum", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "manual_annotations", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
dev.off()

pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "pseudobulk_captureArea_PCA_manual_annotations_4.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "brnum", ncomponents = 4, point_size = 4, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "manual_annotations", ncomponents = 4, point_size = 4, label_format = c("%s %02i", " (%i%%)"),
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
vars <- getVarianceExplained(spe_pseudo, variables=c("brnum", "manual_annotations","sample_id","age","sex"))
head(vars)
# brnum manual_annotations sample_id         age        sex
# ENSG00000241860 17.483652          15.207832  29.35013 0.026739457 2.10585712
# ENSG00000237491 10.728778          14.438564  17.34343 0.009682524 2.39479338
# ENSG00000228794  2.373829          22.392428  13.33559 0.562146602 0.06963826
# ENSG00000225880  2.767235           5.334951  14.32160 0.942949369 0.16608513
# ENSG00000230368 13.919185          12.152200  18.71998 0.295618439 0.25657438
# ENSG00000223764  2.884623          27.920931  10.99554 0.048410489 0.34098590

pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "plot_explanatory_vars_captureArea.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save file
save(spe_pseudo, file = here::here("processed-data", "08_pseudobulk", "manual_annotations", "spe_pseudo_captureArea_manual_annotations.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()