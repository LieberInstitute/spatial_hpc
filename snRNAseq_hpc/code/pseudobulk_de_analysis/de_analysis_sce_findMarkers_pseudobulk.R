#cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(spatialLIBD)
    library(here)
    library(edgeR)
    library(scuttle)
    library(scater)
    library(scran)
    library(dplyr)
    library(PCAtools)
    library(gridExtra)
    library(ggforce)
    library(pheatmap)
    library(scater)
    library(scran)
})

# Load SCE
load(file=here::here('snRNAseq_hpc','processed-data','sce',
                     'sce_final.rda'))

##load palettes
load(file=here::here('plots','snRNAseq_palettes.rda'))
palette = sn.fine.palette

names(palette) <- levels(sce$fine.cell.class)


##cut amy, thal nuclei to avoid obscuring within-HPC marker genes
sce<-sce[,!sce$fine.cell.class %in% c('Amy','Thal','GABA.PENK')]
sce$fine.cell.class<-droplevels(sce$fine.cell.class)
sce$superfine.cell.class<-droplevels(sce$superfine.cell.class)

###remove genes with no counts in any cell prior to DE testing
no_expr <- which(rowSums(counts(sce)) == 0)
sce <- sce[-no_expr, ]
dim(sce)

#########findMarkers()##############
marks.mid<-findMarkers(sce,groups=sce$mid.cell.class,pval.type='all',
                       direction='up')

marks.fine<-findMarkers(sce,groups=sce$fine.cell.class,pval.type='all',
                       direction='up')

marks.superfine<-findMarkers(sce,groups=sce$superfine.cell.class,pval.type='all',
                       direction='up')

####save outputs
save(marks.mid,file=here::here('snRNAseq_hpc','processed-data',
                               'de_analysis','marks_mid.rda'))

save(marks.fine,file=here::here('snRNAseq_hpc','processed-data',
                               'de_analysis','marks_fine.rda'))

save(marks.superfine,file=here::here('snRNAseq_hpc','processed-data',
                               'de_analysis','marks_superfine.rda'))

############pseudobulk analysis###########
## Pseudobulk by spatial superfine.cell.class
sce_pseudo <- aggregateAcrossCells(
    sce,
    DataFrame(
        superfine.cell.class = colData(sce)$superfine.cell.class,
        sample_ID = colData(sce)$sample_ID
    ))
dim(sce_pseudo)
#[1] 35255  1198

##remove pseudobulked samples with very low ncells (here this means # spots)
sce_pseudo <- sce_pseudo[, sce_pseudo$ncells >= 10]
dim(sce_pseudo)
#[1] 35255   827

###remove irrelevant colData cols
colData(sce_pseudo)<-colData(sce_pseudo)[,c('brnum','round','sort','sex','age',
                                            'superfine.cell.class','fine.cell.class','mid.cell.class','broad.cell.class',
                                            'sample_ID','ncells')]

#find a good expression cutoff using edgeR::filterByExpr
rowData(sce_pseudo)$high_expr_group_sample_ID <- filterByExpr(sce_pseudo, group = sce_pseudo$sample_ID)
rowData(sce_pseudo)$high_expr_group_superfine.cell.class <- filterByExpr(sce_pseudo, group = sce_pseudo$superfine.cell.class)

summary(rowData(sce_pseudo)$high_expr_group_sample_ID)
# Mode   FALSE    TRUE
# logical   14349   20906


summary(rowData(sce_pseudo)$high_expr_group_superfine.cell.class)
# Mode   FALSE    TRUE
# logical   11271   23984

with(rowData(sce_pseudo), table(high_expr_group_sample_ID, high_expr_group_superfine.cell.class))
# high_expr_group_sample_ID FALSE  TRUE
# FALSE 11271  3078
# TRUE      0 20906

## Now filter
sce_pseudo <- sce_pseudo[rowData(sce_pseudo)$high_expr_group_sample_ID, ]
dim(sce_pseudo)
#20906   827

x <- edgeR::cpm(edgeR::calcNormFactors(sce_pseudo), log = TRUE, prior.count = 1)
#
## Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(sce_pseudo)))
#
## Fix the column names. DGEList will have samples name as Sample1 Sample2 etc
dimnames(x) <- dimnames(sce_pseudo)
#
## Store the log normalized counts on the SingleCellExperiment object
logcounts(sce_pseudo) <- x


##quick QC
sce_pseudo <- scuttle::addPerCellQC(
    sce_pseudo,
    subsets = list(Mito = which(seqnames(sce_pseudo) == "chrM"))
)
sce_pseudo$det_out<-as.logical(isOutlier(sce_pseudo$detected,type='lower',nmads=3))
sce_pseudo<-sce_pseudo[,sce_pseudo$det_out==F]
dim(sce_pseudo)
# [1] 20906   796
rm(x)

#run PCA
set.seed(12141)
#runPCA(sce_pseudo)
pca <- prcomp(t(assays(sce_pseudo)$logcounts))

message(Sys.time(), " % of variance explained for the top 20 PCs:")
metadata(sce_pseudo)
metadata(sce_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
metadata(sce_pseudo)


pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(sce_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]


# Plot PCA
pdf(file = here::here("snRNAseq_hpc","plots", "de_analysis", "pseudobulk_PCA_sce_final.pdf"), width = 14, height = 14)
plotPCA(sce_pseudo, colour_by = "fine.cell.class", ncomponents = 5, point_size = 1, point_alpha=1,label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)+scale_color_manual(values=palette)+labs(color='fine.cell.class')
plotPCA(sce_pseudo, colour_by = "broad.cell.class", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)+scale_color_manual(values=sn.broad.palette)+labs(color='broad.cell.class')
plotPCA(sce_pseudo, colour_by = "sex", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "age", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "subsets_Mito_percent", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "sum", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "detected", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
dev.off()

##############DE ANALYSIS################
###superfine.cell.class
sce_pseudo$age_scaled<-scales::rescale(sce_pseudo$age,to=c(0,1))
sce_pseudo$superfine.cell.class<-factor(make.names(sce_pseudo$superfine.cell.class))
mod<-registration_model(
    sce_pseudo,
    covars = c('round','sex','age_scaled'),
    var_registration = "superfine.cell.class"
)

cors<-registration_block_cor(
    sce_pseudo,
    mod,
    var_sample_id = "sample_ID"
)

reg<-registration_stats_enrichment(
    sce_pseudo,
    block_cor=cors,
    covars = c('round','sex','age_scaled'),
    var_registration = "superfine.cell.class",
    var_sample_id = "sample_ID",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

rega<-registration_stats_anova(
    sce_pseudo,
    block_cor=cors,
    covars = c('round','sex','age_scaled'),
    var_registration = "superfine.cell.class",
    var_sample_id = "sample_ID",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

regp<-registration_stats_pairwise(
    sce_pseudo,
    block_cor=cors,
    registration_model=mod,
    var_registration = "superfine.cell.class",
    var_sample_id = "sample_ID",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

##save files
save(rega,file=here::here('snRNAseq_hpc','processed-data','de_analysis','sce_anova_stats_superfine.rda'))
save(regp,file=here::here('snRNAseq_hpc','processed-data','de_analysis','sce_pairwise_stats_superfine.rda'))
save(reg,file=here::here('snRNAseq_hpc','processed-data','de_analysis','sce_enrichment_stats_superfine.rda'))

###fine.cell.class
sce_pseudo$age_scaled<-scales::rescale(sce_pseudo$age,to=c(0,1))
sce_pseudo$fine.cell.class<-factor(make.names(sce_pseudo$fine.cell.class))
mod<-registration_model(
    sce_pseudo,
    covars = c('round','sex','age_scaled'),
    var_registration = "fine.cell.class"
)

cors<-registration_block_cor(
    sce_pseudo,
    mod,
    var_sample_id = "sample_ID"
)

reg<-registration_stats_enrichment(
    sce_pseudo,
    block_cor=cors,
    covars = c('round','sex','age_scaled'),
    var_registration = "fine.cell.class",
    var_sample_id = "sample_ID",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

rega<-registration_stats_anova(
    sce_pseudo,
    block_cor=cors,
    covars = c('round','sex','age_scaled'),
    var_registration = "fine.cell.class",
    var_sample_id = "sample_ID",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

regp<-registration_stats_pairwise(
    sce_pseudo,
    block_cor=cors,
    registration_model=mod,
    var_registration = "fine.cell.class",
    var_sample_id = "sample_ID",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

##save files
save(rega,file=here::here('snRNAseq_hpc','processed-data','de_analysis','sce_anova_stats_fine.rda'))
save(regp,file=here::here('snRNAseq_hpc','processed-data','de_analysis','sce_pairwise_stats_fine.rda'))
save(reg,file=here::here('snRNAseq_hpc','processed-data','de_analysis','sce_enrichment_stats_fine.rda'))

###mid.cell.class
sce_pseudo$age_scaled<-scales::rescale(sce_pseudo$age,to=c(0,1))
sce_pseudo$mid.cell.class<-factor(make.names(sce_pseudo$mid.cell.class))
mod<-registration_model(
    sce_pseudo,
    covars = c('round','sex','age_scaled'),
    var_registration = "mid.cell.class"
)

cors<-registration_block_cor(
    sce_pseudo,
    mod,
    var_sample_id = "sample_ID"
)

reg<-registration_stats_enrichment(
    sce_pseudo,
    block_cor=cors,
    covars = c('round','sex','age_scaled'),
    var_registration = "mid.cell.class",
    var_sample_id = "sample_ID",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

rega<-registration_stats_anova(
    sce_pseudo,
    block_cor=cors,
    covars = c('round','sex','age_scaled'),
    var_registration = "mid.cell.class",
    var_sample_id = "sample_ID",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

regp<-registration_stats_pairwise(
    sce_pseudo,
    block_cor=cors,
    registration_model=mod,
    var_registration = "mid.cell.class",
    var_sample_id = "sample_ID",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

##save files
save(rega,file=here::here('snRNAseq_hpc','processed-data','de_analysis','sce_anova_stats_mid.rda'))
save(regp,file=here::here('snRNAseq_hpc','processed-data','de_analysis','sce_pairwise_stats_mid.rda'))
save(reg,file=here::here('snRNAseq_hpc','processed-data','de_analysis','sce_enrichment_stats_mid.rda'))


#######correlation with sestan and SRT data supp figure plots###########
sce_stats<-reg[,c(1:53)]
rownames(sce_stats)<-reg$gene

##load sestan stats (enrichment model only)
load(file=here::here('snRNAseq_hpc','processed-data',
                         'de_analysis','sce_enrichment_stats_sestan.rda'))
sestan_stats<-reg[,c(1:69)]
rownames(sestan_stats)<-reg$gene

##set up data for correlation and heatmaps
colnames(sce_stats)<-gsub(pattern='t_stat_',
                          replacement='',
                          colnames(sce_stats))

colnames(sestan_stats)<-gsub(pattern='t_stat_',
                             replacement='',
                             colnames(sestan_stats))

sestan_stats<-sestan_stats[,order]

i<-intersect(rownames(sce_stats),rownames(sestan_stats))
sestan_stats<-sestan_stats[i,]
sce_stats<-sce_stats[i,]

pdf(file=here::here('plots','figures','supp_figures','supp_fig_heatmap_franjic.pdf'),h=11,w=11)
pheatmap(cor(sestan_stats,sce_stats),cluster_rows=F,
         cluster_cols=F,color=viridis::viridis(100))
dev.off()

##load SRT stats
load(file=here::here('processed-data','08_pseudobulk',
                     'PRECAST','visiumHE_DE_stats_domain.rda'))

spe_stats<-stats$enrichment[,c(1:16)]

colnames(spe_stats)<-gsub(pattern='t_stat_',replacement='',colnames(spe_stats))

i<-intersect(rownames(sce_stats),rownames(spe_stats))
sce_stats<-sce_stats[i,]
spe_stats<-spe_stats[i,]


pdf(file=here::here('plots','figures','supp_figures','supp_fig_heatmap_spatialReg.pdf'),h=5,w=11)
pheatmap(cor(spe_stats,sce_stats),cluster_rows=F,
         cluster_cols=F,color=viridis::viridis(100))
dev.off()

