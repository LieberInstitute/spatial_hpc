################################
# spatial_HPC project
# Pseudobulk SPE by capture area
# Anthony Ramnauth, Dec 15 2022
################################

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

table(spe$sample_id,spe$PRECAST_k15)
#                   1    2    3    4    5    6    7    8    9   10   11   12
#  V10B01-085_A1    1 2998  518    0    1    3    0    2  102    0    0   14
#  V10B01-085_B1  491  585  748  416    1    1  125  178   69  356    8  184
#  V10B01-085_C1    1  210 1414    1   24    3    0  101  329   94  356  837
#  V10B01-085_D1   69   88  876  269    0    4  654   49  107  272  140  679
#  V10B01-086_A1  461 1126 1196    8    1    8    0  226  204  195  276  806
#  V10B01-086_B1  193  177 1585  383    3    1   53  291  226  459   94  347
#  V10B01-086_C1  327  755  434    4    0    4   18  642   79  182   16   64
#  V10B01-086_D1   56   44 1097  172    0    1  818  488  232  267   48  111
#  V11A20-297_A1 1305  299  949   77  196  379   52   52  251  558  125   55
#  V11A20-297_B1  375    9  631   13 1121   77    0  371  292  117  309 1011
#  V11A20-297_C1  133  143  154   63    0    0   59 1311  150   72  120 1221
#  V11A20-297_D1    9   80  241 1121    2    0  280  573  156   54   53  722
#  V11L05-333_A1  565  533  477   93    0    0  195   75  243  122  112 2418
#  V11L05-333_B1   14   30  448  758  188  565  269  775  359  557   71  625
#  V11L05-333_C1 1380  599  918   87    1    3   21  866  203  218   62  198
#  V11L05-333_D1    0    4  931  505  870  589   63   28  335  849   61  295
#  V11L05-335_A1    1    7  862  417  360  267  348  302  312 1113   25  131
#  V11L05-335_B1   12    7 1239   29  692  404   22   96  629 1287   36  107
#  V11L05-335_C1    0   10 1556    5 1146  439    2  103  318  967   18   84
#  V11L05-335_D1    7    2  665  850  596  181  109  223   95  746   59  226
#  V11L05-336_A1    8    4  649  479  684  424  457   76  298  948   67  274
#  V11L05-336_B1  629  207  555  103  582   24  125  559  563  293   39  457
#  V11L05-336_C1    0    0  918   20 1231  351    5  416  258  672   51  168
#  V11L05-336_D1    0   59 1104   15 1226  225    6  172  427  808   32  230
#  V11U08-081_A1   91  158 1836   18  310  618    1  433  317  456   80   75
#  V11U08-081_B1   10    7 1392   21  586  421    2  274  221  427   31   75
#  V11U08-081_C1    1    2  787  458    1    4 1584   29   97  531   14   46
#  V11U08-081_D1  525  896  685   65    1    8  445  688  113  507   37  338
#  V11U08-084_A1  158  678  993   12    3  258    0 1293  411  414   52  267
#  V11U08-084_B1   40   13 1434  687    0    6  149  125  261 1285   39  102
#  V11U08-084_C1    3  162 1487  555  168  292  398   26  206  983   36  337
#  V11U08-084_D1    2    1  948  835   74    0  233   43  107  325   30  157
#
#                  13   14   15
#  V10B01-085_A1    9   27    0
#  V10B01-085_B1   57   62    0
#  V10B01-085_C1   13  396   12
#  V10B01-085_D1    3  128   19
#  V10B01-086_A1   70  109    1
#  V10B01-086_B1   41  122    0
#  V10B01-086_C1   27   60    0
#  V10B01-086_D1   24   67    6
#  V11A20-297_A1   33  126    0
#  V11A20-297_B1   36   70  233
#  V11A20-297_C1   57   37    2
#  V11A20-297_D1   99   68   11
#  V11L05-333_A1   21  135    0
#  V11L05-333_B1  134  197    1
#  V11L05-333_C1   40  141    0
#  V11L05-333_D1  107  202  122
#  V11L05-335_A1  391  126    0
#  V11L05-335_B1   90  107    7
#  V11L05-335_C1  103  211    2
#  V11L05-335_D1  524  183   35
#  V11L05-336_A1   65  163    6
#  V11L05-336_B1   57   99    1
#  V11L05-336_C1  255   75   42
#  V11L05-336_D1  100  124  156
#  V11U08-081_A1  113  115    3
#  V11U08-081_B1  126  172    0
#  V11U08-081_C1   17   82    0
#  V11U08-081_D1   22  117    0
#  V11U08-084_A1   80  197   17
#  V11U08-084_B1  314  119   22
#  V11U08-084_C1   32  220   84
#  V11U08-084_D1  520   94  922

# Rmove CP clusters & NA cluster
speb = spe[, which(spe$PRECAST_k15 != "9")]
speb = speb[, which(speb$PRECAST_k15 != "15")]
speb = speb[, which(speb$PRECAST_k15 != "<NA>")]

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
# 30359   288

##
pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "ncells_captureArea_wo_9-15-NA_Fncells50.pdf"), width = 14, height = 14)
hist(spe_pseudo$ncells, breaks = 200)
boxplot(ncells ~ spe_pseudo$PRECAST_k15, data = colData(spe_pseudo))
dev.off()

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$PRECAST_k15)

summary(rowData(spe_pseudo)$high_expr_group_sample_id)
#   Mode   FALSE    TRUE
#logical   14619   15740

summary(rowData(spe_pseudo)$high_expr_group_cluster)
#   Mode   FALSE    TRUE
#logical   18005   12354

with(rowData(spe_pseudo), table(high_expr_group_sample_id, high_expr_group_cluster))
#high_expr_group_sample_id FALSE  TRUE
#                    FALSE 14619     0
#                    TRUE   3386 12354

## Now filter
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_cluster, ]
dim(spe_pseudo)
#12354   288

# Store the log normalized counts on the spe object
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)

# Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo)))

# Fix the column names. DGEList will have samples name as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)

# Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo) <- x
dim(spe_pseudo)
#12354   288

rm(x)

#run PCA
pca <- prcomp(t(assays(spe_pseudo)$logcounts))

message(Sys.time(), " % of variance explained for the top 20 PCs:")
metadata(spe_pseudo)
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
metadata(spe_pseudo)
# [1] 18.30 12.60  3.31  2.00  1.87  1.78  1.61  1.40  1.39  1.33  1.31  1.28
#[13]  1.23  1.22  1.19  1.16  1.13  1.10  1.08  1.04

pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 18.300 12.600  3.310  2.000  1.870  1.780  1.610  1.400  1.390  1.330
#[11]  1.310  1.280  1.230  1.220  1.190  1.160  1.130  1.100  1.080  1.040
#[21]  1.010  0.999  0.913  0.898  0.870  0.853  0.842  0.809  0.763  0.738
#[31]  0.727  0.714  0.705  0.679  0.645  0.626  0.614  0.573  0.556  0.543
#[41]  0.520  0.511  0.504  0.477  0.461  0.457  0.445  0.434  0.420  0.408

# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "pseudobulk_captureArea_PCA_wo_9-15-NA_Fncells50.pdf"), width = 14, height = 14)
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

pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "pseudobulk_captureArea_PCA_2_wo_9-15-NA_Fncells50.pdf"), width = 14, height = 14)
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

pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "pseudobulk_captureArea_PCA_4_wo_9-15-NA_Fncells50.pdf"), width = 14, height = 14)
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
#                    brnum PRECAST_k15 sample_id        age       sex
#ENSG00000237491  6.060682    21.05548 10.356976 0.21319981 0.3708976
#ENSG00000228794  1.269831    28.25788  6.567467 0.01876500 0.2295285
#ENSG00000187634  2.801131    44.03837  9.033060 0.16295210 0.9765789
#ENSG00000188976  1.414074    26.78634  5.244415 0.85674209 0.2479809
#ENSG00000187961  1.464660    16.29784  6.406836 0.07088862 0.3025624
#ENSG00000188290 19.347874    42.97575 26.453132 2.82726398 8.0550343

pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "variance_captureArea_wo_9-15-NA_Fncells50.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save file
save(spe_pseudo, file = here::here("processed-data", "08_pseudobulk", "PRECAST", "spe_pseudo_captureArea_wo_9-15-NA_Fncells50.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
