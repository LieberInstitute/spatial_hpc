###project PRECAST clusters to VSPG samples
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(here)
library(RcppML)
library(edgeR)
library(ggspavis)

###get rewritten plotVisium()) script
source(file=here::here('code','NMF','plotVisium_rewrite.R'))

##load palettes
load(file=here::here('plots','spatial_palette_final.rda'))

load(file=here::here('processed-data','05_preprocess_batchCorrection',
                     'spe_norm_final.rda'))
spe_norm<-spe
load(file=here::here('processed-data','06_clustering',
     'PRECAST','spe_precast_HE_domain.rda'))
spe2<-spe_norm[,spe_norm$brnum %in% levels(spe_norm$brnum)[c(11:12)]]
spe_init<-spe

##aggregate domains
spe<-SingleCellExperiment(assays=list(counts=logcounts(spe_init)),
                          colData=colData(spe_init))

##quickly remove Thal and Amy from the data
#Thal
spe$thal<-ifelse(spe$sample_id=='Br8325_V11A20-297_B1' & spe$domain %in% c('SUB','SUB.RHP'),T,F)
spe<-spe[,spe$thal==F]
#Amy
spe$amy<-ifelse(spe$sample_id %in% c('Br6432_V10B01-085_C1','Br6432_V10B01-086_B1') & spe$domain %in% c('RHP'),T,F)
spe$amy<-ifelse(spe$sample_id =='Br6423_V10B01-085_A1',T,spe$amy)
spe<-spe[,spe$amy==F]
dim(spe)
#[1] 31483 146523

spe_pseudo <- aggregateAcrossCells(
    spe,
    DataFrame(
        domain = colData(spe)$cluster))

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$domain)

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
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_cluster, ]
dim(spe_pseudo)
#21481    18


counts<-counts(spe_pseudo)

# Compute column sums
col_sums <- colSums(counts)

# Rescale each column
rescaled_counts <- t(apply(counts, 1, function(row) row / col_sums))

# Check if columns now sum to 1
print(colSums(rescaled_counts))

##now project the reformatted counts
set.seed(1029)
i<-intersect(rownames(spe2),rownames(rescaled_counts))
loadings<-rescaled_counts
loadings<-loadings[rownames(loadings) %in% i,]
spe2<-spe2[rownames(spe2) %in% i,]
loadings<-loadings[match(rownames(spe2),rownames(loadings)),]

proj<-project(loadings,logcounts(spe2),L1=0.1)
proj<-t(proj)
colnames(proj)<-levels(spe_init$cluster)

max_column <- function(mat) {
    max_col_names <- colnames(mat)[apply(mat, 1, which.max)]
    return(max_col_names)
}

# Adding the new column to the matrix
# Convert matrix to data frame first
proj_df <- as.data.frame(proj)
proj_df$max_col <- max_column(proj)

# Checking the result
head(proj_df)
spe2$domain<-proj_df$max_col
spe2$domain<-factor(spe2$domain,levels=levels(spe_init$cluster))
levels(spe2$domain)[2]<-"CA2.4"
levels(spe2$domain)[3]<-"CA2.4"
levels(spe2$domain)[3]<-"CA1"
levels(spe2$domain)[4]<-"CA1"



brains <- unique(spe_norm$brnum)


speb <- spe2[, (colData(spe2)$brnum == brains[[2]])]
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

plotVisium(
    speb,
    spots = TRUE,
    fill = 'domain',
    #highlight = "broad.domain",
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    palette = spatial.palette,
    image=F)

levels(spe2$domain)[11]<-'SLM.SGZ'

data<-colData(spe_init)[,c('brnum','sample_id','domain')]
data2<-colData(spe2)[,c('brnum','sample_id','domain')]
data3<-rbind(data,data2)
save(data3,file=here::here('processed-dataHE+VSPG_domain.rda')


###make some plots
pdf(file=here::here('plots','figures','supp_figures','projected_spatial_domains.pdf'),h=7,w=7)
brains <- unique(spe_norm$brnum)


speb <- spe2[, (colData(spe2)$brnum == brains[[11]])]
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

print(plotVisium(
    speb,
    spots = TRUE,
    fill = 'domain',
    highlight='domain',
    values=spatial.palette,
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    palette = spatial.palette,
    image=F))
brains <- unique(spe_norm$brnum)

speb <- spe2[, (colData(spe2)$brnum == brains[[12]])]
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

print(plotVisium(
    speb,
    spots = TRUE,
    fill = 'domain',
    highlight='domain',
    values=spatial.palette,
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    palette = spatial.palette,
    image=F))
dev.off()
