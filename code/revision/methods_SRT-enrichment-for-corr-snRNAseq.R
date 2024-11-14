library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scater)
library(edgeR)
library(scuttle)
library(spatialLIBD)
set.seed(123)

load(here::here('processed-data','06_clustering','PRECAST','spe_precast_HE_domain.rda'))
table(spe$VSPG)


spe<-SingleCellExperiment(assays=list(counts=counts(spe)),colData=colData(spe),rowData=rowData(spe))
spe$amy<-ifelse(spe$sample_id %in%
                  c('Br6423_V10B01-085_C1','Br6432_V10B01-086_B1') &
                  spe$domain %in% c('RHP'),T,F)
spe$amy<-ifelse(spe$sample_id =='Br6423_V10B01-085_A1',T,spe$amy)
spe$thal<-ifelse(spe$sample_id=='Br8325_V11A20-297_B1' & spe$domain %in% c('SUB','SUB.RHP'),T,F)

spe$cluster = as.character(spe$cluster)
spe$cluster = ifelse(spe$amy==TRUE, "Amygdala",spe$cluster)
spe$cluster = ifelse(spe$thal==TRUE, "Thalamus",spe$cluster)

spe_pseudo <- aggregateAcrossCells(
  spe,
  DataFrame(
    cluster = colData(spe)$cluster,
    sample_id = spe$sample_id
  ))
dim(spe_pseudo)
#31483   524 -- domain
#31483   567 -- cluster (plus thal/amy)
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 10]
dim(spe_pseudo)
#31483   422 -- domain, filtered to >=20
#31483   488 -- cluster (plus thal/amy), filtered to >=10

colData(spe_pseudo)<-colData(spe_pseudo)[,c(50,51,48,54,55,56,24,31,32,33)]

rowData(spe_pseudo)$high_expr_group_domain <- filterByExpr(spe_pseudo, group = spe_pseudo$domain)
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_domain, ]
dim(spe_pseudo)
#13455   422 -- domain
#13510   488 -- cluster (plus thal/amy), filtered to >=10

x <- cpm(calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)
identical(rownames(x), rownames(spe_pseudo))
dimnames(x) <- dimnames(spe_pseudo)
logcounts(spe_pseudo) <- x


spe_pseudo$age_scaled<-scales::rescale(spe_pseudo$age,to=c(0,1))
spe_pseudo$slide<-factor(spe_pseudo$slide)
spe_pseudo$cluster<-factor(spe_pseudo$cluster)

mod<-registration_model(
  spe_pseudo,
  covars = c('sex','age_scaled','slide'),
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
  covars = c('sex','age_scaled','slide'),
  var_registration = "cluster",
  var_sample_id = "sample_id",
  gene_ensembl = 'gene_id',
  gene_name = 'gene_name'
)
write.csv(reg, "processed-data/revision/spe_enrichment_stats_cluster-k18-plus-amy-thal.csv", row.names=T)
