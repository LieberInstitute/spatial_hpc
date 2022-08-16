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
})


# Load SPE
load(file = here::here("processed-data", "06_Clustering", "spe_modify.Rdata"))

## Pseudo-bulk for BayesSpace k = 11 results
spe_pseudo <- aggregateAcrossCells(
  spe,
  DataFrame(
    BayesSpace = colData(spe)$BayesSpace_harmony_k11_nrep10000,
    brnum = colData(spe)$brnum
  ))

spe_pseudo$BayesSpace <- factor(spe_pseudo$BayesSpace)

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr <- filterByExpr(spe_pseudo)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$BayesSpace)
rowData(spe_pseudo)$high_expr_group_brain <- filterByExpr(spe_pseudo, group = spe_pseudo$brnum)

summary(rowData(spe_pseudo)$high_expr)
summary(rowData(spe_pseudo)$high_expr_group_cluster)
summary(rowData(spe_pseudo)$high_expr_group_brain)

with(rowData(spe_pseudo), table(high_expr_group_brain, high_expr_group_cluster))

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

message(Sys.time(), " % of variance explained for the top 4 PCs:")
metadata(spe_pseudo)
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(4)])
metadata(spe_pseudo)
# $PCA_var_explained
# [1] 23.400 12.900  3.220  2.010  1.900  1.700  1.450  1.410  1.310  1.140
# [11]  1.110  1.070  1.060  1.020  0.997  0.953  0.941  0.934  0.919  0.900

pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 23.400 12.900  3.220  2.010  1.900  1.700  1.450  1.410  1.310  1.140
# [11]  1.110  1.070  1.060  1.020  0.997  0.953  0.941  0.934  0.919  0.900
# [21]  0.897  0.880  0.869  0.848  0.832  0.819  0.807  0.801  0.769  0.754
# [31]  0.723  0.702  0.689  0.668  0.629  0.617  0.590  0.573  0.552  0.526
# [41]  0.518  0.504  0.493  0.464  0.441  0.423  0.413  0.392  0.375  0.367

# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "pseudobulk_PCAbrain _BayesSpace11.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "brnum", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "BayesSpace", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
dev.off()

## Compute some reduced dims
set.seed(20220423)
spe_pseudo <- scater::runMDS(spe_pseudo, ncomponents = 20)
spe_pseudo <- scater::runPCA(spe_pseudo, name = "runPCA")

####plot explanatory variables ####

#uses linear regression model
vars <- getVarianceExplained(spe_pseudo, variables=c("brnum", "BayesSpace","sample_id","age","sex"))
head(vars)
#                     brnum BayesSpace sample_id         age        sex
# ENSG00000237491  7.154467   14.73066 12.568606 0.007851715 0.23265035
# ENSG00000228794  4.499692   11.83859 10.864590 0.441771049 1.40248742
# ENSG00000188976  2.396955   12.44802  9.619700 0.054827628 0.99269187
# ENSG00000187961  3.337632   11.33855  8.020533 0.476712973 0.25642626
# ENSG00000188290 12.423113   20.15318 17.852229 0.343189579 2.70226293
# ENSG00000187608 12.513031   18.93802 20.764674 0.028097713 0.08933338

pdf(file = here::here("plots","08_pseudobulk", "plot_explanatory_vars_BayesSpace11.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save RDS file
save(spe_pseudo, file = here::here("processed-data", "08_pseudobulk", "spe_pseudo_brain_BayesSpace11.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()