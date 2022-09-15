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

## Pseudo-bulk for mbkmeans k = 17 results
sce <- as(spe, "SingleCellExperiment")
spe_pseudo <- aggregateAcrossCells(
  sce,
  id=DataFrame(
    manual_annotations = colData(sce)$ManualAnnotation,
    brnum = colData(sce)$brnum)
)

spe_pseudo$manual_annotations <- factor(spe_pseudo$manual_annotations)
dim(spe_pseudo)
# [1] 30359    88

spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 200]
dim(spe_pseudo)
#[1] 30359    81

dim(spe_pseudo)
# [1] 30359    80
##
pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "ncells_brain_wo_CP-THAL-CTX_Fncells200.pdf"), width = 14, height = 14)
hist(spe_pseudo$ncells, breaks = 200)
boxplot(ncells ~ spe_pseudo$manual_annotations, data = colData(spe_pseudo))
dev.off()

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_br <- filterByExpr(spe_pseudo, group = spe_pseudo$brnum.1)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$manual_annotations)

summary(rowData(spe_pseudo)$high_expr_group_cluster)
# Mode   FALSE    TRUE 
# logical   15096   15263 

# Mode   FALSE    TRUE 
# logical   14623   15736 

# Mode   FALSE    TRUE 
# logical   14618   15741 

summary(rowData(spe_pseudo)$high_expr_group_br)
# Mode   FALSE    TRUE 
# logical   15096   15263 

# Mode   FALSE    TRUE 
# logical   14858   15501 

# Mode   FALSE    TRUE 
# logical   14846   15513 

with(rowData(spe_pseudo), table(high_expr_group_br, high_expr_group_cluster))
#                    high_expr_group_cluster
# high_expr_group_br FALSE  TRUE
#             FALSE 15096     0
#             TRUE      0 15263

#                    high_expr_group_cluster
# high_expr_group_br FALSE  TRUE
#             FALSE 14623   235
#             TRUE      0 15501

#                    high_expr_group_cluster
# high_expr_group_br FALSE  TRUE
#             FALSE 14618   228
#             TRUE      0 15513

## Now filter
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_cluster, ]
dim(spe_pseudo)
# [1] 15263    88

# [1] 15736    81

# [1] 15741    80

# Store the log normalized counts on the spe object
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)

# Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo)))

# Fix the column names. DGEList will have samples name as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)

# Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo) <- x
dim(spe_pseudo)
# [1] 15263    88
# [1] 15736    81
# [1] 15741    80

rm(x)

#run PCA
pca <- prcomp(t(assays(spe_pseudo)$logcounts))

message(Sys.time(), " % of variance explained for the top 20 PCs:")
metadata(spe_pseudo)
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
metadata(spe_pseudo)
# $PCA_var_explained
# [1] 23.900 15.600  7.740  5.740  4.900  3.980  2.900  2.650  2.480  2.310
# [11]  2.030  1.990  1.410  1.240  1.160  1.140  0.977  0.889  0.795  0.777

# [1] 23.50 10.00  7.06  5.15  4.30  3.87  3.50  2.66  2.41  2.04  1.84  1.67
# [13]  1.59  1.43  1.37  1.24  1.21  1.07  1.03  1.02

# [1] 24.50 10.50  6.79  4.69  4.19  3.73  2.87  2.56  2.27  2.04  1.78  1.71
# [13]  1.53  1.47  1.34  1.29  1.14  1.11  1.09  1.00

pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 23.900 15.600  7.740  5.740  4.900  3.980  2.900  2.650  2.480  2.310
# [11]  2.030  1.990  1.410  1.240  1.160  1.140  0.977  0.889  0.795  0.777
# [21]  0.703  0.676  0.618  0.615  0.548  0.529  0.502  0.485  0.447  0.419
# [31]  0.403  0.385  0.365  0.353  0.335  0.315  0.306  0.295  0.289  0.267
# [41]  0.260  0.252  0.247  0.245  0.237  0.227  0.224  0.212  0.207  0.204

# [1] 23.500 10.000  7.060  5.150  4.300  3.870  3.500  2.660  2.410  2.040
# [11]  1.840  1.670  1.590  1.430  1.370  1.240  1.210  1.070  1.030  1.020
# [21]  0.944  0.880  0.844  0.800  0.757  0.706  0.693  0.659  0.626  0.611
# [31]  0.594  0.573  0.527  0.518  0.505  0.496  0.494  0.470  0.450  0.434
# [41]  0.420  0.409  0.404  0.382  0.367  0.342  0.340  0.334  0.322  0.319

# [1] 24.500 10.500  6.790  4.690  4.190  3.730  2.870  2.560  2.270  2.040
# [11]  1.780  1.710  1.530  1.470  1.340  1.290  1.140  1.110  1.090  1.000
# [21]  0.936  0.901  0.852  0.805  0.752  0.742  0.704  0.671  0.654  0.639
# [31]  0.611  0.563  0.552  0.540  0.529  0.527  0.502  0.479  0.465  0.448
# [41]  0.436  0.431  0.409  0.389  0.365  0.362  0.357  0.343  0.337  0.326

# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "pseudobulk_brain_PCA_wo_CP-THAL-CTX_Fncells200.pdf"), width = 14, height = 14)
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

pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "pseudobulk_brain_PCA_2_wo_CP-THAL-CTX_Fncells200.pdf"), width = 14, height = 14)
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

pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "pseudobulk_brain_PCA_4_wo_CP-THAL-CTX_Fncells200.pdf"), width = 14, height = 14)
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
#all Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
#                         You're computing too large a percentage of total singular values, use a standard svd instead.

####plot explanatory variables ####

#uses linear regression model
vars <- getVarianceExplained(spe_pseudo, variables=c("manual_annotations","brnum","age","sex"))
head(vars)

#                 manual_annotations     brnum        age       sex
# ENSG00000241860          12.452996 34.328676 0.36817394 1.7633193
# ENSG00000237491          11.967989 17.032831 0.14546844 7.1087299
# ENSG00000228794          39.529039  9.857021 0.06645813 1.7446803
# ENSG00000230368           9.612402 38.112190 0.08512377 1.0925235
# ENSG00000223764          44.422727  7.253375 0.82552140 1.7428269
# ENSG00000187634          44.381149  8.474463 0.13871862 0.1331302

#                 manual_annotations     brnum          age       sex
# ENSG00000241860           15.70961 32.410105 0.4653593121 0.1347495
# ENSG00000237491           12.38252 17.106779 2.5433996117 3.9205503
# ENSG00000228794           46.87644  8.914358 0.0733701824 2.1434076
# ENSG00000230368           13.46065 46.004120 0.3186133447 0.1447717
# ENSG00000223764           52.89826 11.163175 1.5724392717 1.2119325
# ENSG00000187634           56.70172 11.890608 0.0006763469 2.1110293

#                  manual_annotations     brnum         age         sex
# ENSG00000241860           18.82599 30.973560 0.179637495 0.508927543
# ENSG00000237491           13.58481 17.145401 2.197709013 3.468695359
# ENSG00000228794           49.58015  9.394223 0.119317750 2.434865411
# ENSG00000230368           12.87758 42.739615 0.099654396 0.009953036
# ENSG00000223764           51.53949 10.034915 1.278760461 0.936305356
# ENSG00000187634           57.17409 11.873145 0.002424467 2.211696922

pdf(file = here::here("plots","08_pseudobulk", "manual_annotations", "variance_brain_wo_CP-THAL-CTX_Fncells200.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save file
save(spe_pseudo, file = here::here("processed-data", "08_pseudobulk", "manual_annotations", "spe_pseudo_brain_wo_CP-THAL-CTX_Fncells200.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()