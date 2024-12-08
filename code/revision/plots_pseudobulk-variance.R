library(SpatialExperiment)
library(spatialLIBD)
#library(here)
library(edgeR)
#library(scuttle)
library(scater)
library(scran)
library(dplyr)
library(PCAtools)
#library(sessioninfo)
#library(gridExtra)
#library(ggforce)
#library(pheatmap)
#BiocManager::install("LieberInstitute/jaffelab")

#code to generate pseudobulk and all dims pulled from here:
#https://github.com/LieberInstitute/spatial_hpc/blob/main/code/08_pseudobulk/PRECAST/create_pseudobulk_spe_captureArea.R
#the pca plots in that code don't generate the plot needed so the plot code itself (also the variance for the metadata variables) are new


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


#pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "histogram_boxplot_domain.pdf"), width = 14, height = 14)
#hist(spe_pseudo$ncells, breaks = 200)
#boxplot(ncells ~ spe_pseudo$domain, data = colData(spe_pseudo))
#dev.off()

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
spe_pseudo <- runPCA(spe_pseudo)
spe_pseudo$sex = as.factor(spe_pseudo$sex)
getVarianceExplained(spe_pseudo, variables=c("domain","slide","sample_id",
                                             "brnum","age","sex",
                                             "sum","detected","subsets_Mito_percent"))
expl.var = c("domain","slide","sample_id",
             "brnum","age","sex",
             "sum","detected","subsets_Mito_percent")
out.mtx = matrix(NA, nrow=nrow(spe_pseudo), ncol=length(expl.var))
rownames(out.mtx) <- rownames(spe_pseudo)
colnames(out.mtx) <- expl.var

for (i in expl.var) {
  if(class(colData(spe_pseudo)[,i])=="factor") {
    colData(spe_pseudo)[,i] <- droplevels(colData(spe_pseudo)[,i])
  }
  out.mtx[,i] = getVarianceExplained(spe_pseudo, variables=i)
}

out.df = tibble::rownames_to_column(as.data.frame(out.mtx), var="gene") %>% 
  tidyr::pivot_longer(cols=expl.var, names_to="mdata", values_to="variance") %>%
  mutate(mdata2= factor(mdata, levels=expl.var, labels=c("spatial domain", "slide", "capture area",
                                                         "donor", "age", "sex", 
                                                         "library size", "genes detected","mito fraction")))

ggplot(out.df, aes(x=variance, color=mdata2))+
  geom_density(key_glyph = draw_key_)+scale_x_log10(limits=c(.001, 100), labels=c(.001,.01,.1,1,10,100),
                               breaks=c(.001,.01,.1,1,10,100))+
  scale_color_brewer(palette="Set1")#+theme_bw()

p1 <- ggplot(out.df, aes(x=variance, color=mdata2))+
  geom_density(size=1, show_guide=FALSE)+
  stat_density(geom="line",position="identity", size=1)+
  scale_x_log10(limits=c(.001, 100), labels=c(.001,.01,.1,1,10,100),
                breaks=c(.001,.01,.1,1,10,100))+
  scale_color_brewer("", palette="Set1")+
  theme(text=element_text(size=12))

pdf(file="plots/revision/supp_pseudobulk-variance.pdf", height=4, width=6)
p1
dev.off()
