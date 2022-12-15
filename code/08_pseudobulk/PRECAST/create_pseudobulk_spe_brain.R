###############################
# spatial_HPC project
# Pseudobulk SPE by brain donor
# Anthony Ramnauth, Dec 15 2022
###############################

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

# Remove CP clusters & NA cluster
speb = spe[, which(spe$PRECAST_k15 != "9")]
speb = speb[, which(speb$PRECAST_k15 != "15")]
speb = speb[, which(speb$PRECAST_k15 != "<NA>")]


## Pseudo-bulk for PRECAST_k15 results
sce <- as(speb, "SingleCellExperiment")
spe_pseudo <- aggregateAcrossCells(
  sce,
  id=DataFrame(
    PRECAST_k15 = colData(sce)$PRECAST_k15,
    brnum = colData(sce)$brnum)
)

spe_pseudo$PRECAST_k15 <- factor(spe_pseudo$PRECAST_k15)

dim(spe_pseudo)
# [1] 30359   117

##
pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "ncells_brain_wo_9-15-NA.pdf"), width = 14, height = 14)
hist(spe_pseudo$ncells, breaks = 200)
boxplot(ncells ~ spe_pseudo$PRECAST_k15, data = colData(spe_pseudo))
dev.off()

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_br <- filterByExpr(spe_pseudo, group = spe_pseudo$brnum.1)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$PRECAST_k15)

summary(rowData(spe_pseudo)$high_expr_group_cluster)
#   Mode   FALSE    TRUE
#logical   16001   14358

summary(rowData(spe_pseudo)$high_expr_group_br)
#   Mode   FALSE    TRUE
#logical   16497   13862

with(rowData(spe_pseudo), table(high_expr_group_br, high_expr_group_cluster))
#                  high_expr_group_cluster
#high_expr_group_br FALSE  TRUE
#             FALSE 16001   496
#             TRUE      0 13862

## Now filter
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_cluster, ]
dim(spe_pseudo)
# 14358   117

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
# [1] 32.90 11.60  3.83  3.05  2.93  2.90  2.75  2.66  2.32  2.26  2.02  1.84
#[13]  1.72  1.55  1.49  1.30  1.27  1.20  1.08  1.07

pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 32.900 11.600  3.830  3.050  2.930  2.900  2.750  2.660  2.320  2.260
#[11]  2.020  1.840  1.720  1.550  1.490  1.300  1.270  1.200  1.080  1.070
#[21]  0.981  0.883  0.841  0.821  0.768  0.707  0.681  0.653  0.614  0.580
#[31]  0.545  0.509  0.481  0.472  0.440  0.403  0.393  0.380  0.369  0.341
#[41]  0.290  0.288  0.281  0.252  0.233  0.221  0.215  0.193  0.179  0.177


# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "pseudobulk_brain_PCA_wo_9-15-NA.pdf"), width = 14, height = 14)
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

pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "pseudobulk_brain_PCA_2_wo_9-15-NA.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "brnum", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "PRECAST_k15", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
dev.off()

pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "pseudobulk_brain_PCA_4_wo_9-15-NA.pdf"), width = 14, height = 14)
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
vars <- getVarianceExplained(spe_pseudo, variables=c("PRECAST_k15","brnum","age","sex"))
head(vars)

#                PRECAST_k15     brnum        age        sex
#ENSG00000241860   21.789916 18.459040 2.26690106 0.07987281
#ENSG00000237491    9.539473  9.898642 1.42454612 0.68208565
#ENSG00000228794   11.506314 12.830054 3.95437636 1.83417852
#ENSG00000230368   21.785085 28.953256 1.42884451 0.31827228
#ENSG00000223764   55.004883  4.634190 0.26167964 1.28067677
#ENSG00000187634   56.391851  3.226901 0.06354897 1.09778703

pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "variance_brain_wo_9-15-NA.pdf"))
plotExplanatoryVariables(vars)
dev.off()

# save file
save(spe_pseudo, file = here::here("processed-data", "08_pseudobulk", "PRECAST", "spe_pseudo_brain_k15_wo_9-15-NA.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()