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
  library(ggforce)
})

# Load SPE
load(file = here::here("processed-data", "06_clustering", "PRECAST", "spe_modify_PRECAST_k16.Rdata"))

# Maddy left the spots that failed QC, discard the failed spots QC'ed by sample_id
spe <- spe[, colData(spe)$discard_auto_id == FALSE]
dim(spe)

table(spe$sample_id,spe$PRECAST_k16)
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

# Rmove CP clusters & NA cluster
#speb = spe[, which(spe$PRECAST_k16 != "9")]
#speb = speb[, which(speb$PRECAST_k16 != "15")]
speb = spe[, !is.na(spe$cluster)]


## Pseudo-bulk for PRECAST_k16
spe_pseudo <- aggregateAcrossCells(
  spe,
  DataFrame(
    cluster = colData(spe)$cluster,
    sample_id = spe$sample_id
  ))

spe_pseudo$cluster <- factor(spe_pseudo$cluster)
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 50]
colData(spe_pseudo)<-colData(spe_pseudo)[,c(1,21,24,31:33,57,81:84)]

dim(spe_pseudo)
# 30432   409

##
pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "histogram_boxplot_precast16.pdf"), width = 14, height = 14)
hist(spe_pseudo$ncells, breaks = 200)
boxplot(ncells ~ spe_pseudo$cluster, data = colData(spe_pseudo))
dev.off()

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$cluster)

summary(rowData(spe_pseudo)$high_expr_group_sample_id)
#Mode   FALSE    TRUE
#logical   16231   14097


summary(rowData(spe_pseudo)$high_expr_group_cluster)
#   Mode   FALSE    TRUE
#logical   17999   12360

with(rowData(spe_pseudo), table(high_expr_group_sample_id, high_expr_group_cluster))
#high_expr_group_sample_id FALSE  TRUE
#                    FALSE 14752     0
#                    TRUE   1479 14097

## Now filter
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_sample_id, ]
dim(spe_pseudo)
#15576   409

#run voom
x<-voom(counts(spe_pseudo), design = 'mod', lib.size = NULL, 
     block = spe_pseudo$batch, correlation = NULL, weights = NULL,
     span = 0.5, plot = FALSE, save.plot = FALSE)
## Store the log normalized counts on the spe object
#x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)
#
## Verify that the gene order hasn't changed
stopifnot(identical(rownames(x$E), rownames(spe_pseudo)))
#
## Fix the column names. DGEList will have samples name as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)
#
## Store the log normalized counts on the SingleCellExperiment object
#logcounts(spe_pseudo) <- x
#dim(spe_pseudo)
##12360   281
#
#rm(x)

#run PCA
set.seed(12141)
#runPCA(spe_pseudo)
pca <- prcomp(t(assays(spe_pseudo)$logcounts))

message(Sys.time(), " % of variance explained for the top 20 PCs:")
metadata(spe_pseudo)
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
metadata(spe_pseudo)
# $PCA_var_explained
# [1] 14.900  8.630  3.840  2.600  2.130  1.940  1.510  1.360  0.975  0.937  0.888  0.798  0.782
# [14]  0.769  0.756  0.729  0.723  0.697  0.671  0.643


pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 14.900  8.630  3.840  2.600  2.130  1.940  1.510  1.360  0.975  0.937  0.888  0.798  0.782
# [14]  0.769  0.756  0.729  0.723  0.697  0.671  0.643  0.620  0.612  0.583  0.574  0.548  0.540
# [27]  0.537  0.531  0.525  0.519  0.507  0.503  0.489  0.485  0.477  0.470  0.466  0.462  0.458
# [40]  0.452  0.446  0.444  0.439  0.430  0.429  0.427  0.422  0.417  0.413  0.412

# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "pseudobulk_PCA_PRECAST16_2comp.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "broad", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "broad2", ncomponents = 2, point_size = 2, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "broad2", ncomponents = 2, point_size = 2, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "subsets_Mito_percent", ncomponents = 2, point_size = 2, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 2, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
dev.off()

pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "pseudobulk_captureArea_PCA_2_wo_9-15-NA_Fncells50.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "brnum", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))
plotPCA(spe_pseudo, colour_by = "detected", ncomponents = 2, point_size = 3, label_format = c("%s %02i", " (%i%%)"),
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
plotPCA(spe_pseudo, colour_by = "PRECAST_k16", ncomponents = 4, point_size = 4, label_format = c("%s %02i", " (%i%%)"),
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

spe_pseudo <- scuttle::addPerCellQC(
    spe_pseudo,
    subsets = list(Mito = which(seqnames(spe_pseudo) == "chrM")),
    BPPARAM = BiocParallel::MulticoreParam(4)
)
spe_pseudo$det_out<-as.logical(isOutlier(spe_pseudo$detected,type='lower',nmads=3))
spe_pseudo<-spe_pseudo[,spe_pseudo$det_out==F]
spe_pseudo$cluster<-droplevels(spe_pseudo$cluster)
spe_pseudo_2$broad_class<-ifelse(spe_pseudo$cluster %in% levels(spe_pseudo2$cluster),'Vasc',
                         ifelse(spe_pseudo$cluster %in% c(1,4,11),'WM',
                           ifelse(spe_pseudo$cluster %in% c(12,17),'Vasc',
                             ifelse(spe_pseudo$cluster %in% c(5,6,9,13),'Neuropil',
                               'Neuron'))))
spe$broad_class<-ifelse(spe$cluster %in% c(12,17),'Vasc',
                               ifelse(spe$cluster %in% c(1,4,11),'WM',
                                      ifelse(spe$cluster %in% c(12,17),'Vasc',
                                             ifelse(spe$cluster %in% c(5,6,9,13),'Neuropil',
                                                    'Neuron'))))



#uses linear regression model
vars <- getVarianceExplained(spe_pseudo, variables=c("PRECAST_k16","sample_id","age","sex",'detected','subsets_Mito_percent'))
head(vars)
#                    brnum PRECAST_k16 sample_id         age         sex
#ENSG00000237491  5.659548    18.37760 11.435987 0.074224263 0.060950740
#ENSG00000228794  1.442367    28.20134  6.077650 0.067342238 0.007492453
#ENSG00000187634  4.437502    41.13168 12.615463 0.378358463 1.966665636
#ENSG00000188976  2.513945    32.18221  8.967045 0.005082043 0.081176015
#ENSG00000187961  1.612837    16.23562  6.775156 0.001670936 0.034610098
#ENSG00000188290 19.135854    43.32915 25.823902 2.695833055 7.628942721

pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "variance_captureArea_wo_9-15-NA_Fncells50.pdf"))
plotExplanatoryVariables(vars)
dev.off()


plotPCA(spe_pseudo2, colour_by = "cluster", ncomponents = 2, point_size = 2, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)+
    geom_mark_ellipse(aes(color = spe$broad,
    label=spe$broad),
    expand = unit(0.5,"mm"),
    label.buffer = unit(-5, 'mm'))

# save file
save(spe_pseudo, file = here::here("processed-data", "08_pseudobulk", "PRECAST", "spe_pseudo_PRECAST_k16.Rdata"))

# Extract the PCA data from spe_pseudo2 into a dataframe
PCAData <- as.data.frame(reducedDim(spe_pseudo2, "PCA"))
# Assuming "broad" is part of colData in spe_pseudo2
PCAData$broad <- colData(spe_pseudo2)$broad
PCAData$cluster<-spe_pseudo2$cluster

# Now create your PCA plot using ggplot directly
pca_plot <- ggplot(PCAData, aes(PC01, PC02)) +
         geom_point(aes(color = cluster),size=) + # Assuming you have a "cluster" column in your PCAData
         scale_color_manual(values=as.vector(palette36.colors(21)), breaks = levels(PCAData$cluster)) +
         labs(x = "PC1", y = "PC2") +
         theme_bw()+theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank())  # Remove minor grid lines

# Add the ellipses
pdf(file=here::here('plots','figures','figure_2','pca_plot.pdf'),h=6,w=8)
pca_plot + geom_mark_ellipse(aes(color = broad, label = broad),
                             expand = unit(0.5,"mm"),
                             label.buffer = unit(-5, 'mm'),
                             show.legend = FALSE,
                             con.type = 'none')
dev.off()



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


mod<-registration_model(
    spe_pseudo,
    covars = c('sex','age_scaled'),
    var_registration = "broad2"
)

cors<-registration_block_cor(
    spe_pseudo,
    mod,
    var_sample_id = "batch"
)

reg<-registration_stats_enrichment(
    spe_pseudo,
    block_cor=cors,
    covars = c('sex','age_scaled'),
    var_registration = "broad2",
    var_sample_id = "batch",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

rega<-registration_stats_anova(
    spe_pseudo,
    block_cor=cors,
    covars = c('sex','age_scaled'),
    var_registration = "broad2",
    var_sample_id = "batch",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

stats<-list(reg,rega)
names(stats)<-c('enrichment','anova')


tstats1 <- stats$enrichment
tstats1 <- tstats1[ , grepl("^t_stat_", colnames(tstats1))]
colnames(tstats1) <- gsub("^t_stat_", "", colnames(tstats1))
rownames(tstats1) <- stats$anova$ensembl

mod2<-registration_model(
    spe_pseudo,
    covars = 'sex',
    var_registration = "broad"
)

cors2<-registration_block_cor(
    spe_pseudo,
    mod2,
    var_sample_id = "sample_id"
)

reg2<-registration_stats_enrichment(
    spe_pseudo,
    block_cor=cors2,
    covars = 'sex',
    var_registration = "broad",
    var_sample_id = "sample_id",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

rega2<-registration_stats_anova(
    spe_pseudo,
    block_cor=cors2,
    covars = 'sex',
    var_registration = "broad",
    var_sample_id = "sample_id",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)


mod<-registration_model(
    sce_pseudo,
    covars = 'sex',
    var_registration = "cellType"
)

cors<-registration_block_cor(
    sce_pseudo,
    mod,
    var_sample_id = "Sample"
)

reg<-registration_stats_enrichment(
    sce_pseudo,
    block_cor=cors,
    covars = 'sex',
    var_registration = "cellType",
    var_sample_id = "Sample",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

rega3<-registration_stats_anova(
    spe_pseudo,
    block_cor=cors3,
    covars = 'sex',
    var_registration = "broad3",
    var_sample_id = "sample_id",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

rega3<-registration_stats_anova(
    spe_pseudo,
    block_cor=cors3,
    covars = 'sex',
    var_registration = "broad3",
    var_sample_id = "sample_id",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)



stats2<-list(reg2,rega2)
names(stats2)<-c('enrichment','anova')


spe_pseudo$broad3<-factor(ifelse(spe_pseudo$cluster %in% c('ML','GCL'),'DG',as.character(spe_pseudo$broad))
)


stats<-registration_wrapper(
    sce,
    var_registration='cellType',
    var_sample_id='Sample',
    covars = 'sex',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    min_ncells = 10,
    pseudobulk_rds_file = NULL
)

sce_pseudo<-registration_pseudobulk(
    sce,
    var_registration='cellType',
    var_sample_id='Sample',
    covars = 'sex',
    min_ncells = 10,
    pseudobulk_rds_file = NULL
)
