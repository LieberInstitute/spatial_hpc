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
load(file = here("processed-data", "06_Clustering", "mbkmeans.Rdata"))

spe$kmeans <- km_res[[13]]$Clusters
table(spe$sample_id,spe$kmeans)
#                  1    2    3    4    5    6    7    8    9   10   11   12
# V10B01-085_A1   97   81  349    0  132   29  482  123   22 1574  211    0
# V10B01-085_B1  219  109  389    0  262   25  390  189   15  675  213    0
# V10B01-085_C1  229  144  332    0  261  189  295  242   54  556  270    0
# V10B01-085_D1  203  133  293    0  187   51  387  231   21  468  306    0
# V10B01-086_A1  345  173  399    0  345   62  404  417   40  761  390    0
# V10B01-086_B1  324  146  436    0  296   63  424  262   38  612  377    0
# V10B01-086_C1  162   82  253    1  132   33  310  301    9  503  208    0
# V10B01-086_D1  183  159  332    0  272   72  375  321   30  583  286    0
# V11A20-297_A1  321  176  557    0  356   51  417  239   24  780  380    0
# V11A20-297_B1  351  129  405    0  282  256  380  438  107  457  360    2
# V11A20-297_C1  447  145   82    0  217   74   89  750   26  137  139    3
# V11A20-297_D1  315   63  403    0  168   61  271  370   35  483  234    2
# V11L05-333_A1  373  158  322    0  186   46  321  784   53  405  333    0
# V11L05-333_B1  306  222  434    0  372   62  281  626   82  429  522    1
# V11L05-333_C1  262  144  475    1  289   78  474  377   41  736  388    3
# V11L05-333_D1  475  257  593    0  479  118  320  222   71  465  552    0
# V11L05-335_A1  456  257  369    0  438   65  385  382   80  385  532    5
# V11L05-335_B1  400  257  500    1  348  117  455  291  211  393  507    4
# V11L05-335_C1  402  227  593    0  406   61  416  296   88  509  606    0
# V11L05-335_D1  318  211  571    0  481   71  450  253   24  500  450    4
# V11L05-336_A1  368  149  590    0  371   47  445  208   44  691  482    2
# V11L05-336_B1  268  146  385    0  279   59  410  329  120  593  367    6
# V11L05-336_C1  339  181  509    0  404   79  454  240   64  602  490    1
# V11L05-336_D1  330  236  452    1  396  170  427  264  116  652  424    0
# V11U08-081_A1  302  183  580    0  317   86  245  223   51  589  561    3
# V11U08-081_B1  248  162  466    0  317   57  223  222   27  389  494    0
# V11U08-081_C1  243  130  443    0  244   34  542  160   12 1011  258    0
# V11U08-081_D1  287  160  454    0  275   44  496  447   18  759  331    1
# V11U08-084_A1  402  263  382    0  425   96  226  654   74  417  525    4
# V11U08-084_B1  530  340  501    1  568   58  340  341   63  365  438    1
# V11U08-084_C1  530  272  482    0  536   99  443  303   30  527  503    1
# V11U08-084_D1  323  192  479    0  545  738  380  131   85  428  335    3
# 
# 13   14   15   16   17
# V10B01-085_A1   83   59   46   83  283
# V10B01-085_B1  195   91   38  166  287
# V10B01-085_C1  402  311   16  158  279
# V10B01-085_D1  309  217   10  191  341
# V10B01-086_A1  400  282   89  202  346
# V10B01-086_B1  248  184   45  146  363
# V10B01-086_C1  159   99   34  127  193
# V10B01-086_D1  225  143   11  153  244
# V11A20-297_A1  207  193  182  183  366
# V11A20-297_B1  533  319   33  266  341
# V11A20-297_C1  601  367   30  201  121
# V11A20-297_D1  409  164   45  181  208
# V11L05-333_A1  752  401   38  428  369
# V11L05-333_B1  462  328  192  233  433
# V11L05-333_C1  280  188   53  339  420
# V11L05-333_D1  299  250  267  175  395
# V11L05-335_A1  238  175  168  229  443
# V11L05-335_B1  230  211  272  184  356
# V11L05-335_C1  176  235  277  178  417
# V11L05-335_D1  298  177  108  220  347
# V11L05-336_A1  234  147  227  178  356
# V11L05-336_B1  346  141   72  231  441
# V11L05-336_C1  267  147  223  155  288
# V11L05-336_D1  284  164  128  206  408
# V11U08-081_A1  146  221  455  227  420
# V11U08-081_B1  140  190  324  140  364
# V11U08-081_C1  103   79   45   98  247
# V11U08-081_D1  327  172   49  178  445
# V11U08-084_A1  408  203  131  278  334
# V11U08-084_B1  252  198   12  225  308
# V11U08-084_C1  295  175  138  259  393
# V11U08-084_D1  220  111    3  131  187

speb = spe[, which(spe$kmeans != "4")]
speb = speb[, which(speb$kmeans != "12")]
speb = speb[, which(speb$kmeans != "6")]
speb = speb[, which(speb$kmeans != "9")]

## Pseudo-bulk for mbkmeans k = 17 results
spe_pseudo <- aggregateAcrossCells(
  speb,
  DataFrame(
    mbkmeans = colData(speb)$kmeans,
    sample_id = speb$sample_id
  ))

spe_pseudo$mbkmeans <- factor(spe_pseudo$mbkmeans)
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 50]

dim(spe_pseudo)
# [1] 30359   502

# [1] 30359   416 (wo clusters 4,12,6,9)

# [1] 30359   401 (wo clusters 4,12,6,9 & filter ncells>50)

##
pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "ncells_captureArea_wo_4-12-6-9_Fncells50.pdf"), width = 14, height = 14)
hist(spe_pseudo$ncells, breaks = 200)
boxplot(ncells ~ spe_pseudo$mbkmeans, data = colData(spe_pseudo))
dev.off()

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$mbkmeans)

summary(rowData(spe_pseudo)$high_expr_group_sample_id)
# Mode   FALSE    TRUE 
# logical   18833   11526 

# Mode   FALSE    TRUE 
# logical   18725   11634 

# Mode   FALSE    TRUE 
# logical   18686   11673 

summary(rowData(spe_pseudo)$high_expr_group_cluster)
# Mode   FALSE    TRUE 
# logical   17130   13229 

# Mode   FALSE    TRUE 
# logical   19550   10809 

# Mode   FALSE    TRUE 
# logical   18903   11456 

with(rowData(spe_pseudo), table(high_expr_group_sample_id, high_expr_group_cluster))
#                           high_expr_group_cluster
# high_expr_group_sample_id FALSE  TRUE
#                     FALSE 17130  1703
#                     TRUE      0 11526

#                           high_expr_group_cluster
# high_expr_group_sample_id FALSE  TRUE
#                     FALSE 18725     0
#                     TRUE    825 10809

#                            high_expr_group_cluster
# high_expr_group_sample_id FALSE  TRUE
#                    FALSE  18686     0
#                    TRUE    217 11456

## Now filter
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_cluster, ]
dim(spe_pseudo)
# [1] 13229   502

# [1] 10809   416 (wo clusters 4,12,6,9)

# [1] 11456   401 (wo clusters 4,12,6,9 & filter ncells>50)
 
# Store the log normalized counts on the spe object
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)

# Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo)))

# Fix the column names. DGEList will have samples name as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)

# Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo) <- x
dim(spe_pseudo)
# [1] 13229   502

# [1] 10809   416 (wo clusters 4,12,6,9)

# [1] 11456   401 
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

# $PCA_var_explained (wo clusters 4,12,6,9)
# [1] 25.500  7.060  5.470  3.750  2.890  1.980  1.650  1.430  1.180  1.060
# [11]  0.906  0.840  0.778  0.717  0.689  0.642  0.607  0.586  0.542  0.519

# $PCA_var_explained 
# [1] 26.200  5.930  4.110  3.300  2.120  1.940  1.680  1.280  1.070  0.966
# [11]  0.934  0.798  0.766  0.704  0.654  0.598  0.596  0.575  0.557  0.533

pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 43.400  8.560  2.630  2.010  1.570  1.220  0.962  0.786  0.639  0.602
# [11]  0.487  0.448  0.443  0.438  0.429  0.425  0.416  0.395  0.389  0.384
# [21]  0.376  0.371  0.367  0.355  0.347  0.342  0.340  0.329  0.327  0.323
# [31]  0.316  0.309  0.307  0.302  0.299  0.290  0.286  0.282  0.276  0.273
# [41]  0.260  0.256  0.251  0.248  0.243  0.239  0.236  0.234  0.232  0.226

# [1] 25.500  7.060  5.470  3.750  2.890  1.980  1.650  1.430  1.180  1.060
# [11]  0.906  0.840  0.778  0.717  0.689  0.642  0.607  0.586  0.542  0.519
# [21]  0.515  0.490  0.475  0.457  0.447  0.441  0.421  0.414  0.400  0.391
# [31]  0.367  0.360  0.347  0.344  0.332  0.327  0.320  0.318  0.312  0.311
# [41]  0.305  0.296  0.293  0.288  0.281  0.277  0.273  0.271  0.266  0.264

# [1] 26.200  5.930  4.110  3.300  2.120  1.940  1.680  1.280  1.070  0.966
# [11]  0.934  0.798  0.766  0.704  0.654  0.598  0.596  0.575  0.557  0.533
# [21]  0.510  0.497  0.474  0.441  0.437  0.424  0.406  0.390  0.389  0.384
# [31]  0.382  0.372  0.365  0.364  0.351  0.347  0.336  0.333  0.332  0.325
# [41]  0.317  0.315  0.314  0.310  0.304  0.302  0.300  0.295  0.292  0.290

# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "pseudobulk_captureArea_PCA_wo_4-12-6-9_Fncells50.pdf"), width = 14, height = 14)
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

pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "pseudobulk_captureArea_PCA_2_wo_4-12-6-9_Fncells50.pdf"), width = 14, height = 14)
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

pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "pseudobulk_captureArea_PCA_4_wo_4-12-6-9_Fncells50.pdf"), width = 14, height = 14)
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
#                 brnum     mbkmeans  sample_id         age         sex
# ENSG00000237491 5.1131169 26.86047  9.033304 0.076173031 0.049237682
# ENSG00000228794 0.9484676 45.58115  3.378215 0.007751482 0.153853642
# ENSG00000187634 5.1216385 32.67619 12.519473 0.299421815 2.886575179
# ENSG00000188976 1.6498779 60.27809  5.290631 0.071596550 0.002581295
# ENSG00000187961 1.4295300 26.71936  5.246689 0.037105746 0.003321983
# ENSG00000188290 6.2041552 60.43831 10.528557 0.038585217 4.018723476

#                   brnum    mbkmeans  sample_id        age        sex
# ENSG00000228794  4.852904 15.81535  7.268836 0.26218853 2.05250725
# ENSG00000188976  3.630736 23.77881  9.440578 0.03093326 0.02385326
# ENSG00000188290 20.594178 35.00307 27.084935 2.01069440 8.85917718
# ENSG00000187608 17.969811 24.21490 27.009651 0.03531768 0.26171027
# ENSG00000188157 14.432904 13.61982 28.284526 0.56086167 0.03784727
# ENSG00000131591  5.520478  9.04122 20.498037 0.12487739 0.84533899

#                     brnum  mbkmeans sample_id        age          sex
# ENSG00000228794  4.397630 20.78298  8.417007 0.76941221 1.7633081843
# ENSG00000187634  4.200752 34.68329 12.792706 0.40717979 1.9356151141
# ENSG00000188976  3.237314 32.20737  8.232627 0.37765592 0.0002202613
# ENSG00000188290 20.837667 36.50443 28.676128 2.03276628 9.2038899767
# ENSG00000187608 21.567532 21.92658 32.500363 0.02107912 0.5261123600
# ENSG00000188157 15.073399 16.24093 30.631928 0.89354011 0.0346158620

pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "variance_captureArea_wo_4-12-6-9_Fncells50.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save file
save(spe_pseudo, file = here::here("processed-data", "08_pseudobulk", "mbkmeans", "spe_pseudo_captureArea_wo_4-12-6-9_Fncells50.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()