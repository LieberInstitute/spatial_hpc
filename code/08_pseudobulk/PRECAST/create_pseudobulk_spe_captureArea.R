#cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/

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
load(file=here::here('processed-data','06_clustering','PRECAST','spe_precast_HE.rda'))
# Set up palette
palette <- c(
  "#008000",  # Level 1 Brighter Forest Green
  "#FF1493",  # Level 2 Bright Hot Pink
  "#FFB6C1",  # Level 3 Medium Pink
  "#800080",  # Level 4 Purple (darker muted purple)
  "#D8BFD8",  # Level 5 Thistle (light muted purple)
  "#FFEC8B",  # Level 6 Light Goldenrod
  "#DAA520",  # Level 7 Goldenrod
  "#FF0000",  # Level 8 Bright Red (less brown, more red)
  "#98FF98",  # Level 9 Mint Green (mintier)
  "#111111",  # Level 10 Lighter Black
  "#dddddd",  # Level 11 Lighter Gray
  "#A9A9A9",  # Level 12 Medium Gray (lighter)
  "#000080",  # Level 13 Lighter Deep Navy
  "#87CEEB",  # Level 14 Sky Blue
  "#4682B4",  # Level 15 Steel Blue
  "#3333DD",  # Level 16 Lighter Dark Blue
  "#F4A460",  # Level 17 Sandy Brown
  "#8B4513"   # Level 18 Saddle Brown
)
names(palette) <- levels(spe$cluster)

# table(spe$sample_id,spe$cluster)
# GCL CA2.4.1 CA2.4.2 CA1.1 CA1.2  SUB SUB.RHP  RHP GABA SL.SR   ML SR.SLM SLM.WM WM.1 WM.2 WM.3 Vascular
# V10B01-085_A1    2       0       1     0     0    0       4 3562   15     2    0      0      0   14    0    2       52
# V10B01-085_B1    1       0       4   112   417  290     786  346   73   153    0    798      7  150   63   27       36
# V10B01-085_C1    0      19      92     0     0    1       8  690  347     7    0    294   1049  811   41   32      345
# V10B01-085_D1    2       2     184   630   330   35     158   75  152   141    0    430    374  672   57   10       95
# V10B01-086_A1    0       0       5     0     3  463     582 1038   83   104  187    494    586  724  314   21       51
# V10B01-086_B1    0       1       7    31   326  128     588  267   95   426  293   1059    281  275   78    4      105
# V10B01-086_C1    4       0       4    16     4  327     285  704   44    87    6    334     16   64  473  199       39
# V10B01-086_D1    1       0       9   773   351   44      67   40   56   679    1    323    160   75  614   45      150
# V11A20-297_A1  188     102     105    20    35 1072     917  197  100   100  228   1029    171   29   19   26       94
# V11A20-297_B1   38     839     595     0     1  290     183    9   69    43   60    407    422  957  281   25      232
# V11A20-297_C1    0       0       1    40    76  129     204   57   41     2    0     28     47  564 2081   98       61
# V11A20-297_D1    0       0       6   232  1169   13      54   56   67   116    0    167     72  675  621   96       67
# V11L05-333_A1    0       0       4   129   164  647     731   87  113    43    0    110    273 2085  507   15       61
# V11L05-333_B1  247     111      99   245   801   22      48    4  207   213  419    810    121  565  646  375       52
# V11L05-333_C1    2       0       1    11    68 1180    1538  117  102    59    6    467     81   87  211  562       56
# V11L05-333_D1  283     450     657    62   477    0      18    2  202   245  645   1055    312  256   25   20      133
# V11L05-335_A1  161     166     384   333   414    2      19    2  133   469  206   1797    209  122   19  100       71
# V11L05-335_B1   44     566     748     5    13   13      18    2   97   179  510   1657    425   51  136  106      166
# V11L05-335_C1  248     518    1155     2     3    0      25    3  221   171  383   1599    327   58   34   66       73
# V11L05-335_D1   96     438     285   134   559    7      27    0  199   906  107    920    119  126  335  144       65
# V11L05-336_A1  218     350     649   366   539   10      25    0  149   120  486   1072    114  247   23   71      100
# V11L05-336_B1   14     378     266    88   124  585     698   93   96    29   21    611    172  394  153  359      112
# V11L05-336_C1  192     615    1084     5    10    0       8    1   80    65  344   1368    166  121  196   62      104
# V11L05-336_D1  126     364    1440     5    17    0      94    5  138    66  233   1416    129  166  121   81      120
# V11U08-081_A1  426     209     536     1    17   95     621   91  131   169  446    971    534   92   47   70      152
# V11U08-081_B1  303     246     992     1    16   10      22   14  180   104  325    962    348   51   70    8      111
# V11U08-081_C1    4       0       4  1574   522    1       2    1   63   628    0    723     14   30   34   15       34
# V11U08-081_D1    5       0       1   411   107  467    1180  150  117    77    3    716     33  174  933   16       53
# V11U08-084_A1  121       0       7     0     2  207    1275  207  190    41  159    533    244  207 1192  283      154
# V11U08-084_B1    5       0       4   112   678   36      60    3  111  1241    2   1976     56   55  105   24       73
# V11U08-084_C1  146     141     375   289   606    2      50  131  211   684  310   1289    247  260   87   15       68
# V11U08-084_D1    0      73      54   224   823    1       5    0  104  1077    0    663     72  100   88   19      121
# V12F14-051_A1  114     381     515   106   258    0       9    1   92   147  166   1290    108  154  250   65       72
# V12F14-051_B1  161     473     726     0     3    0      14   26   96    52  319    866    130   41    5   23      139
# V12F14-051_C1  101      81     236   118   461   43      95    7  153   817  151   1082     96   86  324  113       73
# V12F14-051_D1  242     563     796     0     2    2     168   21  150   132  578   1053    305   70  142   44      169
# 
# Choroid
# V10B01-085_A1       0
# V10B01-085_B1       0
# V10B01-085_C1       2
# V10B01-085_D1       1
# V10B01-086_A1       0
# V10B01-086_B1       0
# V10B01-086_C1       0
# V10B01-086_D1       1
# V11A20-297_A1       0
# V11A20-297_B1     208
# V11A20-297_C1       0
# V11A20-297_D1       1
# V11L05-333_A1       0
# V11L05-333_B1       0
# V11L05-333_C1       0
# V11L05-333_D1      96
# V11L05-335_A1       0
# V11L05-335_B1       1
# V11L05-335_C1       1
# V11L05-335_D1      16
# V11L05-336_A1       0
# V11L05-336_B1       0
# V11L05-336_C1      22
# V11L05-336_D1     137
# V11U08-081_A1       1
# V11U08-081_B1       0
# V11U08-081_C1       0
# V11U08-081_D1       0
# V11U08-084_A1       0
# V11U08-084_B1       0
# V11U08-084_C1      75
# V11U08-084_D1     867
# V12F14-051_A1       0
# V12F14-051_B1       0
# V12F14-051_C1       0
# V12F14-051_D1       1

# Rmove CP clusters & NA cluster
#speb = spe[, which(spe$PRECAST_k16 != "9")]
#speb = speb[, which(speb$PRECAST_k16 != "15")]
#speb = spe[, !is.na(spe$cluster)]


## Pseudo-bulk for PRECAST_k16
spe_pseudo <- aggregateAcrossCells(
  spe,
  DataFrame(
    cluster = colData(spe)$cluster,
    sample_id = spe$sample_id
  ))

spe_pseudo$cluster <- factor(spe_pseudo$cluster)
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 50]
colData(spe_pseudo)<-colData(spe_pseudo)[,c(21,22,24,31:33,48,49,50)]
colData(spe_pseudo)$tissue.type<-factor(
  ifelse(spe_pseudo$cluster %in% levels(spe_pseudo$cluster)[1:9],'Neuron',
         ifelse(spe_pseudo$cluster %in% levels(spe_pseudo$cluster)[10:13],'Neuropil',
                ifelse(spe_pseudo$cluster %in% levels(spe_pseudo$cluster)[14:16],'WM',
                       'Vasc/CSF'))))

dim(spe_pseudo)
# 31483   409

##
pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "histogram_boxplot_precast16.pdf"), width = 14, height = 14)
hist(spe_pseudo$ncells, breaks = 200)
boxplot(ncells ~ spe_pseudo$cluster, data = colData(spe_pseudo))
dev.off()

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$cluster)

summary(rowData(spe_pseudo)$high_expr_group_sample_id)
# Mode   FALSE    TRUE 
# logical   15907   15576


summary(rowData(spe_pseudo)$high_expr_group_cluster)
# Mode   FALSE    TRUE 
# logical   17386   14097

with(rowData(spe_pseudo), table(high_expr_group_sample_id, high_expr_group_cluster))
#high_expr_group_sample_id FALSE  TRUE
#                    FALSE 14752     0
#                    TRUE   1479 14097

## Now filter
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_sample_id, ]
dim(spe_pseudo)
#15576   409

#run voom
# x<-voom(counts(spe_pseudo), design = 'mod', lib.size = NULL, 
#      block = spe_pseudo$batch, correlation = NULL, weights = NULL,
#      span = 0.5, plot = FALSE, save.plot = FALSE)
## Store the log normalized counts on the spe object
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)
#
## Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo)))
#
## Fix the column names. DGEList will have samples name as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)
#
## Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo) <- x


dim(spe_pseudo)
##quick QC
spe_pseudo <- scuttle::addPerCellQC(
  spe_pseudo,
  subsets = list(Mito = which(seqnames(spe_pseudo) == "chrM")),
  BPPARAM = BiocParallel::MulticoreParam(4)
)
spe_pseudo$det_out<-as.logical(isOutlier(spe_pseudo$detected,type='lower',nmads=3))
spe_pseudo<-spe_pseudo[,spe_pseudo$det_out==F]
dim(spe_pseudo)
# [1] 15576   373
rm(x)

#run PCA
set.seed(12141)
#runPCA(spe_pseudo)
pca <- prcomp(t(assays(spe_pseudo)$logcounts))

message(Sys.time(), " % of variance explained for the top 20 PCs:")
metadata(spe_pseudo)
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
metadata(spe_pseudo)
# $PCA_var_explained
# [1] 18.000  5.010  3.880  3.400  2.820  2.260  2.040  1.730  1.220  0.977  0.955  0.869  0.853  0.759  0.665  0.626  0.594
# [18]  0.594  0.578  0.561


pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 18.000  5.010  3.880  3.400  2.820  2.260  2.040  1.730  1.220  0.977  0.955  0.869  0.853  0.759  0.665  0.626  0.594
# [18]  0.594  0.578  0.561  0.553  0.531  0.517  0.509  0.506  0.489  0.485  0.477  0.474  0.471  0.468  0.461  0.457  0.447
# [35]  0.441  0.440  0.436  0.433  0.427  0.418  0.414  0.411  0.407  0.400  0.396  0.392  0.384  0.381  0.379  0.375

# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "pseudobulk_PCA_visiumHE.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "cluster", ncomponents = 4, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)+scale_color_manual(values=palette)
plotPCA(spe_pseudo, colour_by = "tissue.type", ncomponents = 4, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "dateImg", ncomponents = 4, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 4, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 4, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "subsets_Mito_percent", ncomponents = 4, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sum", ncomponents = 4, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "detected", ncomponents = 4, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
dev.off()

vars <- getVarianceExplained(spe_pseudo, variables=c("cluster","sample_id","age_scaled","sex",
                                                     'detected','subsets_Mito_percent','tissue.type','dateImg'))
head(vars)
#                    brnum PRECAST_k16 sample_id         age         sex
#ENSG00000237491  5.659548    18.37760 11.435987 0.074224263 0.060950740
#ENSG00000228794  1.442367    28.20134  6.077650 0.067342238 0.007492453
#ENSG00000187634  4.437502    41.13168 12.615463 0.378358463 1.966665636
#ENSG00000188976  2.513945    32.18221  8.967045 0.005082043 0.081176015
#ENSG00000187961  1.612837    16.23562  6.775156 0.001670936 0.034610098
#ENSG00000188290 19.135854    43.32915 25.823902 2.695833055 7.628942721

pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "variance_explained_visiumHE.pdf"))
plotExplanatoryVariables(vars)
dev.off()

getExplanatoryPCs(spe_pseudo)

# save file
save(spe_pseudo, file = here::here("processed-data", "08_pseudobulk", "PRECAST", "spe_pseudo_PRECAST_k16.Rdata"))

# Extract the PCA data from spe_pseudo into a dataframe
PCAData <- as.data.frame(reducedDim(spe_pseudo, "PCA"))
# Assuming "broad" is part of colData in spe_pseudo
PCAData$tissue.type <- colData(spe_pseudo)$tissue.type
PCAData$cluster<-spe_pseudo$cluster

# Now create your PCA plot using ggplot directly
pca_plot <- ggplot(PCAData, aes(PC01, PC02)) +
         geom_point(aes(color = cluster),size=2) + # Assuming you have a "cluster" column in your PCAData
         scale_color_manual(values=palette, breaks = levels(PCAData$cluster)) +
         labs(x = "PC1", y = "PC2") +
         theme_bw()+theme(panel.grid.major = element_blank(),  # Remove major grid lines
                          panel.grid.minor = element_blank())  # Remove minor grid lines

# Add the ellipses
pdf(file=here::here('plots','figures','figure_2','pca_plot.pdf'),h=6,w=8)
pca_plot + geom_mark_ellipse(aes(color = tissue.type, label = tissue.type),
                             expand = unit(0.5,"mm"),
                             label.buffer = unit(-5, 'mm'),
                             show.legend = FALSE,
                             con.type = 'none')
dev.off()

##############DE ANALYSIS################
spe_pseudo$age_scaled<-scales::rescale(spe_pseudo$age,to=c(0,1))
spe_pseudo$dateImg<-factor(
  gsub(spe_pseudo$dateImg,pattern='-',replacement='_'))
levels(spe_pseudo$tissue.type)[3]<-'Vasc_CSF'
mod<-registration_model(
    spe_pseudo,
    covars = c('sex','age_scaled','dateImg'),
    var_registration = "tissue.type"
)

cors<-registration_block_cor(
    spe_pseudo,
    mod,
    var_sample_id = "sample_id"
)

reg<-registration_stats_enrichment(
    spe_pseudo,
    block_cor=cors,
    covars = c('sex','age_scaled','dateImg'),
    var_registration = "tissue.type",
    var_sample_id = "sample_id",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

rega<-registration_stats_anova(
  spe_pseudo,
  block_cor=cors,
  covars = c('sex','age_scaled','dateImg'),
  var_registration = "tissue.type",
  var_sample_id = "sample_id",
  gene_ensembl = 'gene_id',
  gene_name = 'gene_name'
)

regp<-registration_stats_pairwise(
  spe_pseudo,
  block_cor=cors,
  registration_model=mod,
  var_registration = "tissue.type",
  var_sample_id = "sample_id",
  gene_ensembl = 'gene_id',
  gene_name = 'gene_name'
)

stats<-list(reg,rega,regp)
names(stats)<-c('enrichment','anova','pairwise')
save(stats,file=here::here('processed-data','08_pseudobulk','PRECAST','visiumHE_DE_stats_tissue.rda'))

mod<-registration_model(
  spe_pseudo,
  covars = c('sex','age_scaled','dateImg'),
  var_registration = "cluster"
)

cors<-registration_block_cor(
  spe_pseudo,
  mod,
  var_sample_id = "sample_id"
)

reg<-registration_stats_enrichment(
  spe_pseudo,
  block_cor=cors,
  covars = c('sex','age_scaled','dateImg'),
  var_registration = "cluster",
  var_sample_id = "sample_id",
  gene_ensembl = 'gene_id',
  gene_name = 'gene_name'
)

rega<-registration_stats_anova(
  spe_pseudo,
  block_cor=cors,
  covars = c('sex','age_scaled','dateImg'),
  var_registration = "cluster",
  var_sample_id = "sample_id",
  gene_ensembl = 'gene_id',
  gene_name = 'gene_name'
)

regp<-registration_stats_pairwise(
  spe_pseudo,
  block_cor=cors,
  registration_model=mod,
  var_registration = "cluster",
  var_sample_id = "sample_id",
  gene_ensembl = 'gene_id',
  gene_name = 'gene_name'
)

stats<-list(reg,rega,regp)
names(stats)<-c('enrichment','anova','pairwise')
save(stats,file=here::here('processed-data','08_pseudobulk','PRECAST','visiumHE_DE_stats_cluster.rda'))
save(spe_pseudo,file=here::here('processed-data','08_pseudobulk','PRECAST','spe_pseudo_HE.rda'))
