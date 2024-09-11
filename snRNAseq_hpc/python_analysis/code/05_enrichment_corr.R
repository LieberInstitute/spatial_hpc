library(SingleCellExperiment)
library(edgeR)
library(scuttle)
library(spatialLIBD)

load('snRNAseq_hpc/processed-data/sce/sce_final.rda')
plot(1:10)
sce_pseudo <- aggregateAcrossCells(
  sce,
  DataFrame(
    fine.cell.class = colData(sce)$fine.cell.class,
    sample_ID = colData(sce)$Sample
  ))
dim(sce_pseudo)
#36601 553

sce_pseudo <- sce_pseudo[, sce_pseudo$ncells >= 10]
dim(sce_pseudo)
#36601 448

colData(sce_pseudo)<-colData(sce_pseudo)[,c("Sample",'brnum','round','sort','sex','age',
                                            'fine.cell.class','mid.cell.class',
                                            'ncells')]

rowData(sce_pseudo)$high_expr_group_Sample <- filterByExpr(sce_pseudo, group = sce_pseudo$Sample)
summary(rowData(sce_pseudo)$high_expr_group_Sample)

sce_pseudo <- sce_pseudo[rowData(sce_pseudo)$high_expr_group_Sample, ]
dim(sce_pseudo)

x <- edgeR::cpm(edgeR::calcNormFactors(sce_pseudo), log = TRUE, prior.count = 1)
identical(rownames(x), rownames(sce_pseudo))
dimnames(x) <- dimnames(sce_pseudo)
logcounts(sce_pseudo) <- x

### DE ANALYSIS

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
  var_sample_id = "Sample"
)

reg<-registration_stats_enrichment(
  sce_pseudo,
  block_cor=cors,
  covars = c('round','sex','age_scaled'),
  var_registration = "fine.cell.class",
  var_sample_id = "Sample",
  gene_ensembl = 'gene_id',
  gene_name = 'gene_name'
)

write.csv(reg, "snRNAseq_hpc/python_analysis/processed-data/sce_enrichment_stats_fine.csv", row.names=T)
reg = read.csv("snRNAseq_hpc/python_analysis/processed-data/sce_enrichment_stats_fine.csv", row.names=1)

sce_stats<-reg[,grep("^t_stat",colnames(reg))]
rownames(sce_stats)<-reg$ensembl

save(reg,file=here::here())
load('snRNAseq_hpc/processed-data/de_analysis/sce_enrichment_stats_fine.rda')

##############_ADATA_################
# prep sce from adata if not already done
{
  load("snRNAseq_hpc/python_analysis/processed-data/sce_filtered-matrix.Rdata")
  mdata = colData(sce)
  mdata$old_key = rownames(mdata)
  rownames(mdata) = mdata$key
  
  obs.mtx = read.csv("snRNAseq_hpc/python_analysis/processed-data/adata_qc-strict_scvi-2k-hdg_processed_obs.csv", row.names=3)
  rownames(obs.mtx) = paste(obs.mtx$sample, obs.mtx$barcode, sep="_")
  mdata_filt = mdata[rownames(obs.mtx),]
  
  count.mtx = counts(sce)
  count.mtx = count.mtx[,mdata_filt$old_key]
  colnames(count.mtx) = rownames(mdata_filt)
  identical(rownames(mdata_filt), rownames(obs.mtx))
  
  #var.mtx = read.csv("snRNAseq_hpc/python_analysis/processed-data/adata_qc-strict_scvi-2k-hdg_processed_var.csv", row.names=1)
  #var.mtx$gene_name = rownames(var.mtx)
  #count.mtx = read.csv("snRNAseq_hpc/python_analysis/processed-data/adata_qc-strict_scvi-2k-hdg_processed_counts.csv")
  
  umap.mtx = read.csv("snRNAseq_hpc/python_analysis/processed-data/adata_qc-strict_scvi-2k-hdg_processed_umap.csv", header=F)
  colnames(umap.mtx) = c("UMAP1","UMAP2")
  rownames(umap.mtx) = rownames(obs.mtx)
  
  adata = SingleCellExperiment(list(counts=count.mtx),
                               colData=obs.mtx,
                               rowData=rowData(sce))
  
  reducedDim(adata, "UMAP") = umap.mtx
  counts(adata) <- as(counts(adata), "dgCMatrix")
  save(adata, file="snRNAseq_hpc/python_analysis/processed-data/sce_qc-strict_scvi-2k-hdg_processed.Rdata")
}

load("snRNAseq_hpc/python_analysis/processed-data/sce_qc-strict_scvi-2k-hdg_processed.Rdata")

adata_pseudo <- aggregateAcrossCells(
  adata,
  DataFrame(
    clusters = colData(adata)$scvi_leiden_k10,
    sample_ID = colData(adata)$sample
  ))
dim(adata_pseudo)
#36601 608

adata_pseudo <- adata_pseudo[, adata_pseudo$ncells >= 10]
dim(adata_pseudo)
#36601 527

adata_pseudo <- adata_pseudo[rownames(sce_stats), ]
dim(adata_pseudo)
#21969 527

y <- edgeR::cpm(edgeR::calcNormFactors(adata_pseudo), log = TRUE, prior.count = 1)
identical(rownames(y), rownames(adata_pseudo))
dimnames(y) <- dimnames(adata_pseudo)
logcounts(adata_pseudo) <- y

### DE ANALYSIS

adata_pseudo$age_scaled<-scales::rescale(adata_pseudo$age,to=c(0,1))
adata_pseudo$clusters<-factor(make.names(adata_pseudo$clusters))

mod2 <-registration_model(
  adata_pseudo,
  covars = c('round','sex','age_scaled'),
  var_registration = "clusters"
)

cors2 <-registration_block_cor(
  adata_pseudo,
  mod2,
  var_sample_id = "sample"
)

reg2 <-registration_stats_enrichment(
  adata_pseudo,
  block_cor=cors2,
  covars = c('round','sex','age_scaled'),
  var_registration = "clusters",
  var_sample_id = "sample",
  gene_ensembl = 'ID',
  gene_name = 'Symbol'
)

write.csv(reg2, "snRNAseq_hpc/python_analysis/processed-data/adata_enrichment_stats_scvi-leiden-k10.csv", row.names=T)
reg2 = read.csv("snRNAseq_hpc/python_analysis/processed-data/adata_enrichment_stats_scvi-leiden-k10.csv", row.names=1)

adata_stats<-reg2[,grep("^t_stat",colnames(reg2))]

##############_CORRELATE_ENRICHMENT_RESULTS_################

colnames(sce_stats)<-gsub(pattern='t_stat_',
                          replacement='',
                          colnames(sce_stats))

colnames(adata_stats)<-gsub(pattern='t_stat_X',
                             replacement='c',
                             colnames(adata_stats))

cor.mtx = cor(adata_stats,sce_stats)

pheatmap::pheatmap(cor.mtx[r_ordered,c_ordered],
                   clustering_method = "ward.D",
                   cluster_rows=F, cluster_cols=F,
                   color=viridis::viridis(100))

c_ordered = c("Astro", "Oligo", "OPC", "Choroid", "Ependy", "Vascular", "Micro.Macro.T",
              "Thal","GABA.MGE", "GABA.CGE", "GABA.LAMP5", "GABA.PENK",
              "Amy","HATA", "GC", "CA2.4", "CA1.ProS", "Sub.1",
            "L6.6b", "Sub.2", "L2.3.PrS.Ent", "L5.6", "L2.3.PrS.PaS"
)
setdiff(colnames(cor.mtx), c_ordered)
setdiff(c_ordered, colnames(cor.mtx))
r_ordered = c("c29","c13","c14","c21","c11","c30","c33","c6","c16","c24",
              "c1","c27","c3","c20","c26","c0","c9","c22","c7",
              "c25","c18","c17","c19","c23","c5","c31",
              "c12","c10","c32","c8","c15","c2","c4","c28"
)
setdiff(rownames(cor.mtx), r_ordered)
setdiff(r_ordered, rownames(cor.mtx))

###############
sce_exc = sce[,sce$fine.cell.class %in% c("CA2-4","CA1/ProS","Sub.1","Sub.2","L6/6b",
                                          "L2/3.PrS.Ent","L5/6","L2/3.PrS.PaS")]

library(scater)
plotUMAP(sce_exc, color_by="k_5_louvain")

