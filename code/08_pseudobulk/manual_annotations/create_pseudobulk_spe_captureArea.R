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

# Load SPE
load(file = here::here("processed-data", "06_Clustering", "spe_modify.Rdata"))

dim(spe)
# [1]  30359 135640

table(spe$sample_id,spe$ManualAnnotation)
#                 CA4   CP  CTX  GCL   ML PCL-CA1 PCL-CA3  SGZ   SL  SLM   SO
# V10B01-085_A1    0    0  171    0    0       0       0    0    0    0    0
# V10B01-085_B1    0    0    0    0    0     327       0    0    0  427   22
# V10B01-085_C1    0  268  254    0    0     397       0    0    0 1283  265
# V10B01-085_D1    0   51    0    0    0     500       0    0    0  873  824
# V10B01-086_A1    0    0    0    0    0     890       0    0    0 1411  130
# V10B01-086_B1    0   88    0    0    0    1426       0    0    0  611  929
# V10B01-086_C1    0    0    0    0    0     393       0    0    0  247  280
# V10B01-086_D1    0    0    0    0    0    1486       0    0    0  177  436
# V11A20-297_A1  219    0    0  193  460     817       0   86    0  581  113
# V11A20-297_B1  556  489    0   46   93       0     449   69   85  320  475
# V11A20-297_C1    0    0    0    0    0     335       0    0    0    0  184
# V11A20-297_D1    0    0    0    0    0    1436       0    0    0   33  395
# V11L05-333_A1    0    0    0    0    0     813       0    0    0    0  301
# V11L05-333_B1  287    0    0  235  304    1074       0  104    0  650  723
# V11L05-333_C1    0    0    0    0    0     786       0    0    0    0    0
# V11L05-333_D1  894   88    0  295  689     641     249  168  320  463  282
# V11L05-335_A1  760    0    0  223  277    1006       0  150 1093    3    0
# V11L05-335_B1    0   86    0  249  385     130     743    0    0 1085  602
# V11L05-335_C1 2090   48    0  341  386       0     246  175  795    0  219
# V11L05-335_D1  470   38    0  118  222    1499     246  156  372  201  349
# V11L05-336_A1  820    0    0  278  585     828      76  403  165  807  105
# V11L05-336_B1  585  101    0   43   26       0       0  230    0  850    0
# V11L05-336_C1  863   69    0  344  622       0     653  396  395  584  266
# V11L05-336_D1 1326  185    0  171  238       0     493  271  228 1145  213
# V11U08-081_A1  759   38    0  672  537    1018     131  654   31  328   30
# V11U08-081_B1 1484    0    0  408  529       0      18  417    0  284    0
# V11U08-081_C1    0    0    0    0    0    2075       0    0    0  295   37
# V11U08-081_D1    0    0    0    0    0     800       0    0    0  218  563
# V11U08-084_A1   67    0    0  145  230       0       0   96    0  350  389
# V11U08-084_B1    0    0    0    8  112     983       0    0    0 1059  255
# V11U08-084_C1   45   82    0  216  598    1081     150  154   29  950  485
# V11U08-084_D1    0 1027    0    0    0     638     517    0  118 1196  440
# 
#                 SR  SUB THAL   WM
# V10B01-085_A1    0 3483    0    0
# V10B01-085_B1    0 1769    0  718
# V10B01-085_C1  432    0    0  834
# V10B01-085_D1  467  532    0   98
# V10B01-086_A1  410 1812    0    0
# V10B01-086_B1  324  513    0   73
# V10B01-086_C1  308  828    0  533
# V10B01-086_D1  782    0    0  502
# V11A20-297_A1  674 1125    0  158
# V11A20-297_B1    0    0  304 1773
# V11A20-297_C1    0    0    0 2910
# V11A20-297_D1  118    0    0 1430
# V11L05-333_A1  183  256    0 3416
# V11L05-333_B1  465    0    0 1143
# V11L05-333_C1  645 1795    0 1322
# V11L05-333_D1  570    0    0  254
# V11L05-335_A1 1082    0    0    0
# V11L05-335_B1 1132    0    0  188
# V11L05-335_C1  576    0    0    0
# V11L05-335_D1  587    0    0  218
# V11L05-336_A1  331    0    0  141
# V11L05-336_B1    0 1725    0  633
# V11L05-336_C1    0    0    0  251
# V11L05-336_D1    0  122    0  266
# V11U08-081_A1  411    0    0    0
# V11U08-081_B1  517  106    0    0
# V11U08-081_C1 1219    0    0   23
# V11U08-081_D1  859 1190    0  813
# V11U08-084_A1  147 2122    0 1276
# V11U08-084_B1 1136  873    0  115
# V11U08-084_C1  747  171    0  278
# V11U08-084_D1  238    0    0  117

spe = spe[, which(spe$ManualAnnotation != "CP")]
spe = spe[, which(spe$ManualAnnotation != "THAL")]
spe = spe[, which(spe$ManualAnnotation != "CTX")]

dim(spe)
# [1]  30359 132021

## Pseudo-bulk for manualannotations k = 15 results
spe_pseudo <- aggregateAcrossCells(
  spe,
  DataFrame(
    manual_annotations = colData(spe)$ManualAnnotation,
    sample_id = spe$sample_id
  ))

spe_pseudo$manual_annotations <- factor(spe_pseudo$manual_annotations)
dim(spe_pseudo)
# [1] 30359   232

spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 50]
dim(spe_pseudo)
# [1] 30359   218

##
pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "ncells_captureArea_wo_CP-THAL-CTX_Fncells50.pdf"), width = 14, height = 14)
hist(spe_pseudo$ncells, breaks = 200)
boxplot(ncells ~ spe_pseudo$manual_annotations, data = colData(spe_pseudo))
dev.off()

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$manual_annotations)

summary(rowData(spe_pseudo)$high_expr_group_sample_id)
# Mode   FALSE    TRUE 
# logical   13454   16905 

# Mode   FALSE    TRUE 
# logical   14371   15988 

summary(rowData(spe_pseudo)$high_expr_group_cluster)
# Mode   FALSE    TRUE 
# logical   17275   13084 

# Mode   FALSE    TRUE 
# logical   17094   13265 

with(rowData(spe_pseudo), table(high_expr_group_sample_id, high_expr_group_cluster))
#                           high_expr_group_cluster
# high_expr_group_sample_id FALSE  TRUE
#                     FALSE 13454     0
#                     TRUE   3821 13084

#                           high_expr_group_cluster
# high_expr_group_sample_id FALSE  TRUE
#                     FALSE 14371     0
#                     TRUE   2723 13265

## Now filter
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_cluster, ]
dim(spe_pseudo)
# [1] 13084   232

# [1] 13265   218

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
# [1] 18.400 13.700  4.640  3.190  2.610  2.400  2.320  2.280  2.050  1.860
# [11]  1.640  1.580  1.490  1.320  1.260  1.120  1.070  0.957  0.904  0.877

# $PCA_var_explained
# [1] 20.100  6.840  5.330  3.700  3.100  2.470  2.030  1.770  1.460  1.410
# [11]  1.340  1.220  1.180  1.170  1.100  1.040  0.964  0.921  0.905  0.834

pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 18.400 13.700  4.640  3.190  2.610  2.400  2.320  2.280  2.050  1.860
# [11]  1.640  1.580  1.490  1.320  1.260  1.120  1.070  0.957  0.904  0.877
# [21]  0.823  0.800  0.769  0.740  0.713  0.701  0.671  0.637  0.601  0.587
# [31]  0.578  0.521  0.512  0.502  0.484  0.475  0.463  0.461  0.441  0.431
# [41]  0.430  0.419  0.406  0.395  0.382  0.377  0.375  0.369  0.352  0.350

# [1] 20.100  6.840  5.330  3.700  3.100  2.470  2.030  1.770  1.460  1.410
# [11]  1.340  1.220  1.180  1.170  1.100  1.040  0.964  0.921  0.905  0.834
# [21]  0.823  0.798  0.777  0.762  0.745  0.736  0.709  0.691  0.678  0.655
# [31]  0.645  0.615  0.611  0.607  0.598  0.572  0.567  0.553  0.527  0.524
# [41]  0.502  0.480  0.476  0.471  0.452  0.448  0.434  0.420  0.413  0.403

# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "pseudobulk_captureArea_PCA_wo_CP-THAL-CTX_Fncells50.pdf"), width = 14, height = 14)
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

pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "pseudobulk_captureArea_PCA_2_wo_CP-THAL-CTX_Fncells50.pdf"), width = 14, height = 14)
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

pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "pseudobulk_captureArea_PCA_4_wo_CP-THAL-CTX_Fncells50.pdf"), width = 14, height = 14)
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
#                     brnum manual_annotations sample_id        age        sex
# ENSG00000237491  9.926138           9.813934  17.99382 0.30188917 2.82533724
# ENSG00000228794  4.335046          13.251490  11.54376 0.93045284 1.45203140
# ENSG00000187634  2.520136          30.561013  12.25834 0.15761455 0.53033032
# ENSG00000188976  5.805161          17.775417  18.49064 0.04201076 0.36246970
# ENSG00000187961  4.079720           5.459124  16.63414 0.12265524 0.03935035
# ENSG00000188290 15.469259          27.413786  24.48699 1.46820462 8.00987948

# brnum manual_annotations sample_id         age        sex
# ENSG00000237491 14.243660          10.942856  23.44580 0.783246884  3.7148063
# ENSG00000228794  4.515593          27.337936  12.01316 0.003742368  1.1296760
# ENSG00000187634  5.151541          37.589960  15.51457 0.113766377  1.1046661
# ENSG00000188976  5.408228          37.207475  14.41255 0.202601195  1.0899026
# ENSG00000187961  2.282289           7.413336  21.20220 0.165908867  0.6095876
# ENSG00000188290 22.063032          43.512749  29.05320 1.010382896 10.5379294

pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "variance_captureArea_wo_CP-THAL-CTX_Fncells50.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save file
save(spe_pseudo, file = here::here("processed-data", "08_pseudobulk", "manual_annotations", "spe_pseudo_captureArea_wo_CP-THAL-CTX_Fncells50.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()