discard<-c(259:264)
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
    library(ggspavis)
})

# Load SPE
load(file = here::here("processed-data", "06_clustering", "PRECAST", "spe_modify_PRECAST_k16.Rdata"))

# Maddy left the spots that failed QC, discard the failed spots QC'ed by sample_id
spe <- spe[, colData(spe)$discard_auto_id == FALSE]
dim(spe)

spe$batch<-factor(
  ifelse(spe$brnum %in% levels(spe$brnum)[c(11:12)],'VSPG','H&E'))

table(spe$sample_id,spe$PRECAST_k24)

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

colData(spe_pseudo)<-colData(spe_pseudo)[,c(1,21:24,31:33,218:223)]

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

#ComBat-seq
logcts<-ComBat_seq(
  counts(spe_pseudo),
  batch=spe_pseudo$batch,
  group = spe_pseudo$cluster
)

#run voom
spe_pseudo2<-spe_pseudo
counts(spe_pseudo2)<-logcts
x<-calcNormFactors(spe_pseudo)
x<-voom(x$counts, design = mod,
     block = spe_pseudo$dateImg,
     span = 0.5, plot = T, save.plot = FALSE)
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
#dim(spe_pseudo)
##12360   281
#
#rm(x)
spe_pseudo <- scuttle::addPerCellQC(
    spe_pseudo,
    subsets = list(Mito = which(seqnames(spe_pseudo) == "chrM")),
    BPPARAM = BiocParallel::MulticoreParam(4)
)
spe_pseudo$det_out<-as.logical(isOutlier(spe_pseudo$detected,type='lower',nmads=4))
#discard_index<-rownames(as.data.frame(reducedDim(spe_pseudo,'PCA'))[order(as.data.frame(reducedDim(spe_pseudo,'PCA'))$PC03,decreasing=T),])[1:6]
#spe_pseudo$cp_outliers<-ifelse(colnames(spe_pseudo) %in% discard_index,T,F)
spe_pseudo<-spe_pseudo[,spe_pseudo$det_out==F]

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
plotPCA(spe_pseudo, colour_by = 'TTR', ncomponents = 2, point_size = 2, label_format = c("%s %02i", " (%i%%)"),
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
plotPCA(spe_pseudo, colour_by = "de", ncomponents = 2, point_size = 8, label_format = c("%s %02i", " (%i%%)"),
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


spe_pseudo<-spe_pseudo[,spe_pseudo$cp_outliers==F]
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
PCAData <- as.data.frame(reducedDim(spe_pseudo, "PCA"))
# Assuming "broad" is part of colData in spe_pseudo2
PCAData$broad <- colData(spe_pseudo)$broad.class
PCAData$cluster<-spe_pseudo$cluster

# Now create your PCA plot using ggplot directly
pca_plot <- ggplot(PCAData, aes(PC01, PC02)) +
         geom_point(aes(color = broad),size=2) + # Assuming you have a "cluster" column in your PCAData
         scale_color_manual(values=as.vector(palette36.colors(5)), breaks = levels(PCAData$cluster)) +
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
    covars = c('sex','age_scaled'),
    var_registration = 'cluster',
    var_sample_id = "sample_id",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

rega<-registration_stats_anova(
    spe_pseudo,
    block_cor=cors,
    covars = c('sex','age_scaled'),
    var_registration = "cluster",
    var_sample_id = "sample_id",
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
    covars = c('age_scaled','sex'),
    var_registration = "broad.class"
)

cors2<-registration_block_cor(
    spe_pseudo,
    mod2,
    var_sample_id = "batch2"
)

reg2<-registration_stats_enrichment(
    spe_pseudo,
    block_cor=cors2,
    covars = c('age_scaled','sex'),
    var_registration = "broad.class",
    var_sample_id = "batch2",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

rega2<-registration_stats_anova(
    spe_pseudo,
    block_cor=cors2,
    covars = c('age_scaled','sex'),
    var_registration = "broad.class",
    var_sample_id = "batch2",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)


mod<-registration_model(
    spe_pseudo,
    covars = c('sex','age_scaled'),
    var_registration = "cluster"
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
