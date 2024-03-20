########################################
# spatial_HPC project
# Pseudobulk SPE by capture area keep CP
# Anthony Ramnauth, Feb 15 2022
########################################

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
load(file = here::here("processed-data", "06_clustering", "PRECAST", "spe_modify_PRECAST_k15.Rdata"))

# Maddy left the spots that failed QC, discard the failed spots QC'ed by sample_id
spe <- spe[, colData(spe)$discard_auto_id == FALSE]
dim(spe)

table(spe$sample_id,spe$PRECAST_k15)
#                 1    2    3    4    5    6    7    8    9   10   11   12
#  V10B01-085_A1    1 2973  516    0    1    3    0    0   97    0    0   14
#  V10B01-085_B1  491  585  748  416    1    1  125  164   69  355    8  184
#  V10B01-085_C1    1  210 1398    1   24    3    0   59  325   82  356  835
#  V10B01-085_D1   69   87  862  269    0    4  653   39  105  267  140  679
#  V10B01-086_A1  461 1118 1178    8    1    8    0  223  202  195  276  806
#  V10B01-086_B1  193  177 1585  383    3    1   53  269  225  448   94  347
#  V10B01-086_C1  327  754  388    4    0    4   18  615   76  117   16   62
#  V10B01-086_D1   56   44 1081  172    0    1  818  451  232  243   48  111
#  V11A20-297_A1 1305  299  949   77  196  379   52   32  249  548  125   55
#  V11A20-297_B1  375    9  606   13 1114   73    0  370  292  111  309 1011
#  V11A20-297_C1  127  142  117   49    0    0   58 1244  145   67  120 1219
#  V11A20-297_D1    9   80  241 1111    2    0  280  556  152   54   53  721
#  V11L05-333_A1  501  509  325   85    0    0  195   62  194   87  111 2411
#  V11L05-333_B1   14   30  447  758  188  564  269  771  356  539   71  625
#  V11L05-333_C1 1380  586  903   87    1    3   21  679  196  205   62  197
#  V11L05-333_D1    0    4  931  505  870  587   63   20  331  849   61  295
#  V11L05-335_A1    1    7  862  417  360  267  348  299  312 1113   25  131
#  V11L05-335_B1   12    7 1239   29  692  404   22   87  627 1287   36  107
#  V11L05-335_C1    0   10 1556    5 1146  439    2   73  316  966   18   84
#  V11L05-335_D1    7    2  665  850  596  181  109  212   95  734   59  226
#  V11L05-336_A1    8    4  649  479  684  424  457   38  297  945   67  274
#  V11L05-336_B1  629  207  544  103  582   24  125  472  556  290   39  457
#  V11L05-336_C1    0    0  918   20 1231  351    5  402  258  672   51  168
#  V11L05-336_D1    0   59 1104   15 1226  225    6  148  427  803   32  230
#  V11U08-081_A1   91  158 1836   18  310  618    1  422  315  402   79   75
#  V11U08-081_B1   10    7 1387   21  586  421    2  272  221  411   31   75
#  V11U08-081_C1    1    2  786  458    1    4 1584   27   96  530   14   46
#  V11U08-081_D1  525  896  666   65    1    8  445  683  113  448   36  338
#  V11U08-084_A1  158  678  987   12    3  256    0 1291  410  400   52  267
#  V11U08-084_B1   40   13 1434  687    0    6  149   87  260 1285   39  102
#  V11U08-084_C1    3  162 1487  555  168  292  398   24  206  983   36  337
#  V11U08-084_D1    2    1  948  835   74    0  233   43  107  325   30  157
#
#                  13   14   15
#  V10B01-085_A1    0   27    0
#  V10B01-085_B1   51   62    0
#  V10B01-085_C1    7  396   12
#  V10B01-085_D1    2  128   19
#  V10B01-086_A1   42  109    1
#  V10B01-086_B1   33  122    0
#  V10B01-086_C1   19   60    0
#  V10B01-086_D1   17   67    6
#  V11A20-297_A1   28  126    0
#  V11A20-297_B1   29   70  233
#  V11A20-297_C1   23   37    2
#  V11A20-297_D1   62   68   11
#  V11L05-333_A1   16  127    0
#  V11L05-333_B1  118  197    1
#  V11L05-333_C1   36  141    0
#  V11L05-333_D1   98  202  122
#  V11L05-335_A1  339  126    0
#  V11L05-335_B1   74  107    7
#  V11L05-335_C1   58  211    2
#  V11L05-335_D1  504  183   35
#  V11L05-336_A1   41  163    6
#  V11L05-336_B1   46   99    1
#  V11L05-336_C1  247   75   42
#  V11L05-336_D1   98  124  156
#  V11U08-081_A1   99  115    3
#  V11U08-081_B1  110  172    0
#  V11U08-081_C1   15   82    0
#  V11U08-081_D1   21  117    0
#  V11U08-084_A1   67  197   17
#  V11U08-084_B1  298  119   22
#  V11U08-084_C1   31  220   84
#  V11U08-084_D1  520   94  922

# Rmove NA cluster
speb = spe[, which(spe$PRECAST_k15 != "<NA>")]

## Pseudo-bulk for PRECAST_k15
spe_pseudo <- aggregateAcrossCells(
  speb,
  DataFrame(
    PRECAST_k15 = colData(speb)$PRECAST_k15,
    sample_id = speb$sample_id
  ))

spe_pseudo$PRECAST_k15 <- factor(spe_pseudo$PRECAST_k15)
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 50]

dim(spe_pseudo)
# 30359   318

##
pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "ncells_captureArea_NA_Fncells50.pdf"), width = 14, height = 14)
hist(spe_pseudo$ncells, breaks = 200)
boxplot(ncells ~ spe_pseudo$PRECAST_k15, data = colData(spe_pseudo))
dev.off()

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$PRECAST_k15)

summary(rowData(spe_pseudo)$high_expr_group_sample_id)
#   Mode   FALSE    TRUE
#logical   15645   14714

summary(rowData(spe_pseudo)$high_expr_group_cluster)
#   Mode   FALSE    TRUE
#logical   16452   13907

with(rowData(spe_pseudo), table(high_expr_group_sample_id, high_expr_group_cluster))
#high_expr_group_sample_id FALSE  TRUE
#                    FALSE 15645     0
#                    TRUE   807 13907

## Now filter
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_cluster, ]
dim(spe_pseudo)
#13907   318

# Store the log normalized counts on the spe object
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)

# Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo)))

# Fix the column names. DGEList will have samples name as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)

# Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo) <- x
dim(spe_pseudo)
#13907   318

rm(x)

#run PCA
pca <- prcomp(t(assays(spe_pseudo)$logcounts))

message(Sys.time(), " % of variance explained for the top 20 PCs:")
metadata(spe_pseudo)
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
metadata(spe_pseudo)
# [1] 17.100 10.100  3.860  2.540  2.200  1.880  1.730  1.540  1.190  1.170
#[11]  1.140  1.060  1.050  1.020  1.000  0.979  0.920  0.905  0.893  0.868

pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 17.100 10.100  3.860  2.540  2.200  1.880  1.730  1.540  1.190  1.170
#[11]  1.140  1.060  1.050  1.020  1.000  0.979  0.920  0.905  0.893  0.868
#[21]  0.852  0.811  0.801  0.792  0.765  0.752  0.741  0.719  0.677  0.661
#[31]  0.647  0.628  0.613  0.587  0.577  0.559  0.546  0.542  0.524  0.498
#[41]  0.478  0.476  0.466  0.456  0.449  0.438  0.421  0.415  0.405  0.404

# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "pseudobulk_captureArea_PCA_NA_Fncells50.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "brnum", ncomponents = 12, point_size = 3, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "PRECAST_k15", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 12, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
dev.off()

pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "pseudobulk_captureArea_PCA_2_NA_Fncells50.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "brnum", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))
plotPCA(spe_pseudo, colour_by = "PRECAST_k15", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))
plotPCA(spe_pseudo, colour_by = "ncells", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))
plotPCA(spe_pseudo, colour_by = "pmi", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))
plotPCA(spe_pseudo, colour_by = "experimenterSeq", ncomponents = 2, point_size = 8,
    label_format = c("%s %02i", " (%i%%)"), percentVar = metadata(spe_pseudo)$PCA_var_explained) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))
plotPCA(spe_pseudo, colour_by = "slide", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))
dev.off()

pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "pseudobulk_captureArea_PCA_4_NA_Fncells50.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "brnum", ncomponents = 4, point_size = 4, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "PRECAST_k15", ncomponents = 4, point_size = 4, label_format = c("%s %02i", " (%i%%)"),
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
vars <- getVarianceExplained(spe_pseudo, variables=c("brnum", "PRECAST_k15","sample_id","age","sex"))
head(vars)
#              brnum PRECAST_k15 sample_id          age         sex
#LINC01409  5.380993    17.49252 12.785752 0.0189452061 0.009256790
#LINC01128  1.275064    27.56404  6.131802 0.0253776853 0.002689864
#SAMD11     5.540942    35.93907 13.221336 0.2286953491 2.893733188
#NOC2L      2.392905    31.26615  9.541107 0.0007613315 0.016351754
#KLHL17     1.152672    14.93024  4.314553 0.0036940985 0.003570580
#HES4      15.965175    46.96103 22.924233 1.3370432758 7.049929095

pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "variance_captureArea_NA_Fncells50.pdf"))
plotExplanatoryVariables(vars)
dev.off()



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
