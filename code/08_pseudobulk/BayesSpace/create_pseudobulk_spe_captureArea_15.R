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
dir_rdata <- here::here("processed-data", "08_pseudobulk", "BayesSpace")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

# Create directory for pseudobulked plots
dir_plots <- here::here("plots", "08_pseudobulk", "BayesSpace")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

# Load SPE
load(file = here::here("processed-data", "06_clustering", "BayesSpace", "1st_run",
    "BayesSpace_rerun_k15.Rdata"))

# Remove CP clusters (keep in mind when plotting to correct name titles)
speb = spe[, which(spe$BayesSpace_harmony_k15_nrep10000 != "13")]
speb = speb[, which(speb$BayesSpace_harmony_k15_nrep10000 != "14")]
speb = speb[, which(speb$BayesSpace_harmony_k15_nrep10000 != "15")]

## Pseudo-bulk for BayesSpace k = 15 results
spe_pseudo <- aggregateAcrossCells(
  speb,
  DataFrame(
    BayesSpace = colData(speb)$BayesSpace_harmony_k15_nrep10000,
    sample_id = speb$sample_id
  ))

spe_pseudo$BayesSpace <- factor(spe_pseudo$BayesSpace)
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 50]

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr <- filterByExpr(spe_pseudo)
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$BayesSpace)

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
# [1] 24.500  5.750  4.410  3.430  2.570  2.350  2.040  1.420  1.250  1.130
#[11]  1.070  0.991  0.919  0.865  0.803  0.772  0.762  0.691  0.675  0.653

pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 24.500  5.750  4.410  3.430  2.570  2.350  2.040  1.420  1.250  1.130
#[11]  1.070  0.991  0.919  0.865  0.803  0.772  0.762  0.691  0.675  0.653
#[21]  0.627  0.588  0.560  0.548  0.522  0.513  0.494  0.482  0.471  0.446
#[31]  0.432  0.419  0.412  0.407  0.396  0.391  0.379  0.376  0.365  0.361
#[41]  0.353  0.349  0.347  0.338  0.335  0.328  0.322  0.314  0.312  0.307

# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "BayesSpace",
    "pseudobulk_captureArea_PCA_BayesSpace15_wo131415.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "brnum", ncomponents = 12, point_size = 3, label_format = c("%s %02i", " (%i%%)"),
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

pdf(file = here::here("plots","08_pseudobulk", "BayesSpace",
    "pseudobulk_captureArea_PCA_BayesSpace15_2_wo131415.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "brnum", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "BayesSpace", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
dev.off()

pdf(file = here::here("plots","08_pseudobulk", "BayesSpace",
    "pseudobulk_captureArea_PCA_BayesSpace15_4_wo131415.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "brnum", ncomponents = 4, point_size = 4, label_format = c("%s %02i", " (%i%%)"),
            percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "BayesSpace", ncomponents = 4, point_size = 4, label_format = c("%s %02i", " (%i%%)"),
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
vars <- getVarianceExplained(spe_pseudo, variables=c("brnum", "BayesSpace","sample_id","age","sex"))
head(vars)
#                    brnum BayesSpace sample_id          age         sex
#ENSG00000237491  8.813819   9.045857 17.699872 0.3955747087 0.006978852
#ENSG00000228794  1.968772  26.166520  7.195776 0.0000193499 1.052335777
#ENSG00000187634  3.300818  33.456932 15.948618 0.2261248686 1.945396384
#ENSG00000188976  2.413076  25.097592  9.744639 0.4015336640 0.075846875
#ENSG00000188290 21.972497  33.070365 30.155934 1.9276292880 9.967170735
#ENSG00000187608 20.426963  23.406549 32.708423 0.2249512459 0.338779961

pdf(file = here::here("plots","08_pseudobulk", "BayesSpace",
    "plot_explanatory_vars_captureArea_BayesSpace15_wo131415.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save file
save(spe_pseudo, file = here::here("processed-data", "08_pseudobulk", "BayesSpace", "spe_pseudo_captureArea_BayesSpace15_wo131415.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
