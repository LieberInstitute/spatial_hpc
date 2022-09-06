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

## Pseudo-bulk for mbkmeans k = 17 results
sce <- as(spe, "SingleCellExperiment")
spe_pseudo <- aggregateAcrossCells(
  sce,
  id=DataFrame(
    manual_annotations = colData(sce)$ManualAnnotation,
    brnum = colData(sce)$brnum)
)

spe_pseudo$manual_annotations <- factor(spe_pseudo$manual_annotations)
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 10]

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_br <- filterByExpr(spe_pseudo, group = spe_pseudo$brnum.1)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$manual_annotations)
summary(rowData(spe_pseudo)$high_expr_group_cluster)
summary(rowData(spe_pseudo)$high_expr_group_br)
with(rowData(spe_pseudo), table(high_expr_group_br, high_expr_group_cluster))

## Now filter
dim(spe_pseudo)
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
# [1] 17.600 13.700  6.030  5.010  4.530  4.000  3.340  3.140  2.870  2.300
# [11]  2.020  1.940  1.800  1.720  1.560  1.500  1.320  1.130  1.050  0.999

pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 17.600 13.700  6.030  5.010  4.530  4.000  3.340  3.140  2.870  2.300
# [11]  2.020  1.940  1.800  1.720  1.560  1.500  1.320  1.130  1.050  0.999
# [21]  0.862  0.837  0.782  0.748  0.702  0.658  0.642  0.634  0.609  0.589
# [31]  0.551  0.538  0.497  0.489  0.468  0.451  0.449  0.434  0.422  0.396
# [41]  0.391  0.370  0.356  0.350  0.341  0.330  0.323  0.308  0.300  0.295

# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "pseudobulk_brain_PCA.pdf"), width = 14, height = 14)
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

pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "pseudobulk_brain_PCA_2.pdf"), width = 14, height = 14)
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

pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "pseudobulk_brain_PCA_4.pdf"), width = 14, height = 14)
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
vars <- getVarianceExplained(spe_pseudo, variables=c("manual_annotations","brnum","age","sex"))
head(vars)

# manual_annotations     brnum          age          sex
# ENSG00000241860          17.964874 34.293107 0.1138736530 2.172512e+00
# ENSG00000237491          20.094781 14.724833 0.2492243809 2.376426e+00
# ENSG00000228794          42.808340  8.799720 0.2036091236 4.967157e-01
# ENSG00000225880           9.847936  5.334630 0.9348403531 7.402876e-05
# ENSG00000230368          15.369057 30.459859 0.0002583493 1.493305e-02
# ENSG00000223764          41.124735  8.327073 0.8089869439 1.676746e+00

pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "plot_explanatory_vars_brain.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save file
save(spe_pseudo, file = here::here("processed-data", "08_pseudobulk", "manual_annotations", "spe_pseudo_brain.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()