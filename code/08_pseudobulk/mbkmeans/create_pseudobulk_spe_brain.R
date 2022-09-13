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
table(spe$sample_id,spe$kmeans)
# 1    2    3    4    5    6    7    8    9   10   11   12
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
sce <- as(speb, "SingleCellExperiment")
spe_pseudo <- aggregateAcrossCells(
  sce,
  id=DataFrame(
    mbkmeans = colData(sce)$kmeans,
    brnum = colData(sce)$brnum)
)

spe_pseudo$mbkmeans <- factor(spe_pseudo$mbkmeans)
dim(spe_pseudo)
# [1] 30359   147

# [1] 30359   117  (wo clusters 4,12,6,9)

##
pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "ncells_brain_wo_4-12-6-9.pdf"), width = 14, height = 14)
hist(spe_pseudo$ncells, breaks = 200)
boxplot(ncells ~ spe_pseudo$mbkmeans, data = colData(spe_pseudo))
dev.off()

# spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 10]

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_br <- filterByExpr(spe_pseudo, group = spe_pseudo$brnum.1)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$mbkmeans)

summary(rowData(spe_pseudo)$high_expr_group_cluster)
# Mode   FALSE    TRUE 
# logical   15235   15124 

# Mode   FALSE    TRUE 
# logical   15738   14621 

summary(rowData(spe_pseudo)$high_expr_group_br)
# Mode   FALSE    TRUE 
# logical   16594   13765 

# Mode   FALSE    TRUE 
# logical   16112   14247 

with(rowData(spe_pseudo), table(high_expr_group_br, high_expr_group_cluster))
#                    high_expr_group_cluster
# high_expr_group_br FALSE  TRUE
#             FALSE 15235  1359
#             TRUE      0 13765

#                     high_expr_group_cluster
# high_expr_group_br FALSE  TRUE
#             FALSE 15738   374
#             TRUE      0 14247

## Now filter
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_cluster, ]
dim(spe_pseudo)
# [1] 15124   147

# [1] 14621   117 (wo clusters 4,12,6,9)

# Store the log normalized counts on the spe object
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)

# Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo)))

# Fix the column names. DGEList will have samples name as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)

# Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo) <- x

rm(x)

#run PCA
pca <- prcomp(t(assays(spe_pseudo)$logcounts))

message(Sys.time(), " % of variance explained for the top 20 PCs:")
metadata(spe_pseudo)
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
metadata(spe_pseudo)
# $PCA_var_explained
# [1] 68.100  5.680  2.010  1.490  1.290  1.170  1.050  0.978  0.936  0.896
# [11]  0.853  0.799  0.725  0.674  0.653  0.552  0.500  0.466  0.440  0.408

# $PCA_var_explained  (wo clusters 4,12,6,9)
# [1] 33.300  7.260  5.910  3.970  3.190  2.940  2.100  1.760  1.660  1.580
# [11]  1.480  1.370  1.320  1.290  1.230  1.160  1.110  0.958  0.941  0.809

pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 68.100  5.680  2.010  1.490  1.290  1.170  1.050  0.978  0.936  0.896
# [11]  0.853  0.799  0.725  0.674  0.653  0.552  0.500  0.466  0.440  0.408
# [21]  0.373  0.338  0.324  0.313  0.309  0.308  0.287  0.285  0.280  0.267
# [31]  0.260  0.248  0.241  0.226  0.220  0.216  0.204  0.196  0.190  0.187
# [41]  0.174  0.162  0.156  0.149  0.137  0.119  0.115  0.112  0.105  0.104

# [1] 33.300  7.260  5.910  3.970  3.190  2.940  2.100  1.760  1.660  1.580
# [11]  1.480  1.370  1.320  1.290  1.230  1.160  1.110  0.958  0.941  0.809
# [21]  0.751  0.727  0.688  0.670  0.635  0.586  0.580  0.560  0.549  0.526
# [31]  0.521  0.506  0.499  0.481  0.466  0.460  0.451  0.428  0.419  0.414
# [41]  0.404  0.397  0.387  0.382  0.368  0.360  0.352  0.351  0.340  0.334

# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "pseudobulk_brain_PCA_wo_4-12-6-9.pdf"), width = 14, height = 14)
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

pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "pseudobulk_brain_PCA_2_wo_4-12-6-9.pdf"), width = 14, height = 14)
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

pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "pseudobulk_brain_PCA_4_wo_4-12-6-9.pdf"), width = 14, height = 14)
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
vars <- getVarianceExplained(spe_pseudo, variables=c("mbkmeans","brnum","age","sex"))
head(vars)

# mbkmeans     brnum         age       sex
# ENSG00000241860 47.30329 21.073382 0.005758258 2.5406981
# ENSG00000237491 60.52541  6.826028 0.001834419 0.8945204
# ENSG00000228794 93.27664  1.464072 0.101931629 0.1325683
# ENSG00000230368 33.00999 21.690064 2.809853149 0.6668615
# ENSG00000223764 52.18557  4.135560 0.243995288 1.4884063
# ENSG00000187634 65.15615  7.494925 0.639563646 2.8006216

# mbkmeans     brnum        age       sex
# ENSG00000241860 13.878066 37.943766 0.10462220 3.4603217
# ENSG00000237491 22.645493 19.561571 0.47877187 5.9227151
# ENSG00000228794 46.312154  7.578292 0.02792873 3.9096482
# ENSG00000230368  7.772907 39.896436 2.51481818 1.2135245
# ENSG00000223764 50.617266  1.447651 0.26806261 0.4459086
# ENSG00000187634 60.049032  9.343432 0.78164156 3.8024693

pdf(file = here::here("plots","08_pseudobulk", "mbkmeans", "variance_brain_wo_4-12-6-9.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save file
save(spe_pseudo, file = here::here("processed-data", "08_pseudobulk", "mbkmeans", "spe_pseudo_brain_wo_4-12-6-9.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()