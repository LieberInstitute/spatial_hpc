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
    library(pheatmap)
    library(scater)
    library(scran)
})

# Load SPE
load(file=here::here('processed-data','06_clustering',
                     'PRECAST','spe_precast_HE_domain.rda'))

##load palettes
load(file=here::here('plots','spatial_palette_final.rda'))
palette = spatial.palette

names(palette) <- levels(spe$domain)

## Pseudobulk by spatial domain
spe_pseudo <- aggregateAcrossCells(
    spe,
    DataFrame(
        domain = colData(spe)$domain,
        sample_id = colData(spe)$sample_id
    ))
dim(spe_pseudo)
#[1] 31483   524

##remove pseudobulked samples with very low ncells (here this means # spots)
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 50]
dim(spe_pseudo)
#[1] 31483   372

###remove irrelevant colData cols
colData(spe_pseudo)<-colData(spe_pseudo)[,c(21,22,24,31:33,50,52,53,54)]


pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "histogram_boxplot_domain.pdf"), width = 14, height = 14)
hist(spe_pseudo$ncells, breaks = 200)
boxplot(ncells ~ spe_pseudo$domain, data = colData(spe_pseudo))
dev.off()

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_domain <- filterByExpr(spe_pseudo, group = spe_pseudo$domain)

summary(rowData(spe_pseudo)$high_expr_group_sample_id)
# Mode   FALSE    TRUE
# logical   16079   15404


summary(rowData(spe_pseudo)$high_expr_group_domain)
# Mode   FALSE    TRUE
# logical   17569   13914

with(rowData(spe_pseudo), table(high_expr_group_sample_id, high_expr_group_domain))
#high_expr_group_sample_id FALSE  TRUE
#                    FALSE 16079     0
#                    TRUE   1490 13914

## Now filter
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_sample_id, ]
dim(spe_pseudo)
#15404   372

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


##quick QC
spe_pseudo <- scuttle::addPerCellQC(
    spe_pseudo,
    subsets = list(Mito = which(seqnames(spe_pseudo) == "chrM"))
)
spe_pseudo$det_out<-as.logical(isOutlier(spe_pseudo$detected,type='lower',nmads=2))
spe_pseudo<-spe_pseudo[,spe_pseudo$det_out==F]
dim(spe_pseudo)
# [1] 15404   323
rm(x)

##quickly remove Thal and Amy from the data
#Thal
spe_pseudo$thal<-ifelse(spe_pseudo$sample_id=='Br8325_V11A20-297_B1' & spe_pseudo$domain %in% c('SUB','SUB.RHP'),T,F)
spe_pseudo<-spe_pseudo[,spe_pseudo$thal==F]
#Amy
spe_pseudo$amy<-ifelse(spe_pseudo$sample_id %in%
                           c('Br6423_V10B01-085_C1','Br6432_V10B01-086_B1') &
                           spe_pseudo$domain %in% c('RHP'),T,F)
spe_pseudo$amy<-ifelse(spe_pseudo$sample_id =='Br6423_V10B01-085_A1',T,spe_pseudo$amy)
spe_pseudo<-spe_pseudo[,spe_pseudo$amy==F]
dim(spe_pseudo)
#[1] 15404   317

#run PCA
set.seed(12141)
#runPCA(spe_pseudo)
pca <- prcomp(t(assays(spe_pseudo)$logcounts))

message(Sys.time(), " % of variance explained for the top 20 PCs:")
metadata(spe_pseudo)
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
metadata(spe_pseudo)
#$PCA_var_explained
# [1] 18.200  5.890  3.860  3.140  3.000  2.430  2.080  1.650  1.220  1.100  1.040  0.891  0.851  0.803
# [15]  0.726  0.665  0.650  0.619  0.603  0.587


pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]
# [1] 18.200  5.890  3.860  3.140  3.000  2.430  2.080  1.650  1.220  1.100  1.040  0.891  0.851  0.803
# [15]  0.726  0.665  0.650  0.619  0.603  0.587  0.572  0.568  0.552  0.546  0.536  0.526  0.521  0.512
# [29]  0.510  0.500  0.497  0.486  0.484  0.482  0.474  0.465  0.462  0.458  0.452  0.451  0.443  0.440
# [43]  0.436  0.430  0.427  0.425  0.419  0.415  0.410  0.402

# Plot PCA
pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "pseudobulk_PCA_visiumHE_final.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "domain", ncomponents = 5, point_size = 1, point_alpha=1,label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)+scale_color_manual(values=palette)+labs(color='domain')
plotPCA(spe_pseudo, colour_by = "broad.domain", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)+scale_color_manual(values=spatial.palette2)+labs(color='broad.domain')
plotPCA(spe_pseudo, colour_by = "dateImg", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sex", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "age", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "subsets_Mito_percent", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sum", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "detected", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
dev.off()

# Extract the PCA data from spe_pseudo into a dataframe
PCAData <- as.data.frame(reducedDim(spe_pseudo, "PCA"))
PCAData$broad.domain <- colData(spe_pseudo)$broad.domain
PCAData$domain<-spe_pseudo$domain

# Now create your PCA plot using ggplot directly
pca_plot <- ggplot(PCAData, aes(PC01, PC02)) +
    geom_point(aes(color = domain),size=2) + # Assuming you have a "domain" column in your PCAData
    scale_color_manual(values=spatial.palette, breaks = levels(PCAData$domain)) +
    labs(x = "PC1", y = "PC2") +
    theme_bw()+theme(panel.grid.major = element_blank(),  # Remove major grid lines
                     panel.grid.minor = element_blank())  # Remove minor grid lines

# Add the ellipses
pdf(file=here::here('plots','figures','figure_2','pca_plot_final_check.pdf'),h=6,w=8)
pca_plot + geom_mark_ellipse(aes(color = broad.domain, label = broad.domain),
                             expand = unit(0.5,"mm"),
                             label.buffer = unit(-5, 'mm'),
                             show.legend = FALSE,
                             con.type = 'none')
dev.off()

##############DE ANALYSIS################
spe_pseudo$age_scaled<-scales::rescale(spe_pseudo$age,to=c(0,1))
spe_pseudo$slide<-factor(spe_pseudo$slide)
mod<-registration_model(
    spe_pseudo,
    covars = c('sex','age_scaled','slide'),
    var_registration = "broad.domain"
)

cors<-registration_block_cor(
    spe_pseudo,
    mod,
    var_sample_id = "sample_id"
)

reg<-registration_stats_enrichment(
    spe_pseudo,
    block_cor=cors,
    covars = c('sex','age_scaled','slide'),
    var_registration = "broad.domain",
    var_sample_id = "sample_id",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

rega<-registration_stats_anova(
    spe_pseudo,
    block_cor=cors,
    covars = c('sex','age_scaled','slide'),
    var_registration = "broad.domain",
    var_sample_id = "sample_id",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

regp<-registration_stats_pairwise(
    spe_pseudo,
    block_cor=cors,
    registration_model=mod,
    var_registration = "broad.domain",
    var_sample_id = "sample_id",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

stats<-list(reg,rega,regp)
names(stats)<-c('enrichment','anova','pairwise')
save(stats,file=here::here('processed-data','08_pseudobulk','PRECAST','visiumHE_DE_stats_broad_final.rda'))


mod<-registration_model(
    spe_pseudo,
    covars = c('sex','age_scaled','slide'),
    var_registration = "domain"
)

cors<-registration_block_cor(
    spe_pseudo,
    mod,
    var_sample_id = "sample_id"
)

reg<-registration_stats_enrichment(
    spe_pseudo,
    block_cor=cors,
    covars = c('sex','age_scaled','slide'),
    var_registration = "domain",
    var_sample_id = "sample_id",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

rega<-registration_stats_anova(
    spe_pseudo,
    block_cor=cors,
    covars = c('sex','age_scaled','slide'),
    var_registration = "domain",
    var_sample_id = "sample_id",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

regp<-registration_stats_pairwise(
    spe_pseudo,
    block_cor=cors,
    registration_model=mod,
    var_registration = "domain",
    var_sample_id = "sample_id",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

stats<-list(reg,rega,regp)
names(stats)<-c('enrichment','anova','pairwise')
save(stats,file=here::here('processed-data','08_pseudobulk','PRECAST','visiumHE_DE_stats_domain.rda'))
save(spe_pseudo,file=here::here('processed-data','08_pseudobulk','PRECAST','spe_pseudo_HE.rda'))

###make supplementary figures
vars <- getVarianceExplained(spe_pseudo, variables=c("domain",'broad.domain',"age","sex","brnum","slide","sample_id",
                                                     'dateImg','detected','sum','ncells','subsets_Mito_percent'))
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


covs<-c("domain","age","sex","brnum","slide","sample_id",
        'broad.domain','dateImg','detected','sum','ncells','subsets_Mito_percent')
colData(spe_pseudo)<-colData(spe_pseudo)[,covs]
pheatmap(getExplanatoryPCs(spe_pseudo)[,covs])
