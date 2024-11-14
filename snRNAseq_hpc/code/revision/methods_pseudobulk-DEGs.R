library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scater)
set.seed(123)

load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))
unique(sce$brnum)

load("snRNAseq_hpc/processed-data/pseudobulk/sce_enrichment_stats.rda")
#stats_enrichment: t stats, FDR, logFC
s1 = colnames(stats_enrichment)[grep("t_stat",colnames(stats_enrichment))]
s2 = sapply(strsplit(s1, split="_"), function(x) x[[3]])
setdiff(unique(sce$superfine.cell.class), s2)
#the only ones excluded are: Thal, PENK, Ahi, Cajal. included is HATA
# these were performed on super fine clusters but i won't be able to justify excluding AHi and not HATA without doing something different

group_by(as.data.frame(colData(sce)), Sample, superfine.cell.class) %>% 
  tally() %>% filter(n>=10) #864
group_by(as.data.frame(colData(sce)), brnum, superfine.cell.class) %>% 
  tally() %>% filter(n>=10) #432
group_by(as.data.frame(colData(sce)), round, superfine.cell.class) %>% 
       tally() %>% filter(n>=10) #303

################ DE from scratch ################
sce_pseudo <- aggregateAcrossCells(
  sce,
  DataFrame(
    superfine.cell.class = colData(sce)$superfine.cell.class,
    sample_ID = colData(sce)$Sample
  ))
dim(sce_pseudo)
#36601 1276

sce_pseudo <- sce_pseudo[, sce_pseudo$ncells >= 10]
dim(sce_pseudo)
#36601 864

colData(sce_pseudo)<-colData(sce_pseudo)[,c("Sample",'brnum','round','sort','sex','age',
                                            'k_5_louvain','superfine.cell.class','fine.cell.class',
                                            'ncells')]

rowData(sce_pseudo)$high_expr_group_Sample <- filterByExpr(sce_pseudo, group = sce_pseudo$Sample)
summary(rowData(sce_pseudo)$high_expr_group_Sample)

sce_pseudo <- sce_pseudo[rowData(sce_pseudo)$high_expr_group_Sample, ]
dim(sce_pseudo)
#21104   864

x <- edgeR::cpm(edgeR::calcNormFactors(sce_pseudo), log = TRUE, prior.count = 1)
identical(rownames(x), rownames(sce_pseudo))
dimnames(x) <- dimnames(sce_pseudo)
logcounts(sce_pseudo) <- x

### DE ANALYSIS

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
  var_sample_id = "Sample"
)

reg<-registration_stats_enrichment(
  sce_pseudo,
  block_cor=cors,
  covars = c('round','sex','age_scaled'),
  var_registration = "superfine.cell.class",
  var_sample_id = "Sample",
  gene_ensembl = 'gene_id',
  gene_name = 'gene_name'
)

write.csv(reg, "snRNAseq_hpc/processed-data/revision/sn_enrichment_stats_superfine.csv", row.names=T)
