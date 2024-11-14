library(SingleCellExperiment)
library(edgeR)
library(scuttle)
library(spatialLIBD)
#i guess i didn't set seed >:[

##############_ADATA_################
# prep sce from adata if not already done
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

#reload if already done
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
