library("SingleCellExperiment")
library("jaffelab")
library("scater")
library("scran")
library("here")
library("sessioninfo")

load(here("snRNAseq_hpc","processed-data", "sce", "sce_harmony.rda"), verbose = TRUE)

message("running buildSNNGraph - ", Sys.time())
snn.gr <- buildSNNGraph(sce, k = 20, use.dimred = "HARMONY")
snn.gr2 <- buildSNNGraph(sce, k = 50, use.dimred = "HARMONY")

message("running walktrap - ", Sys.time())
clust <- igraph::cluster_walktrap(snn.gr)$membership
clust50 <- igraph::cluster_walktrap(snn.gr2)$membership

table(clust)
table(clust50)

##add to sce
sce$k_20_label<-clust
sce$k_50_label<-clust50

##make some prelim plots
pdf(here("snRNAseq_hpc","plots","UMAP_k20.pdf"))
plotUMAP(sce,colour_by='k_20_label',text_by='k_20_label')
dev.off()

pdf(here("snRNAseq_hpc","plots","UMAP_k50.pdf"))
plotUMAP(sce,colour_by='k_50_label',text_by='k_50_label')
dev.off()

##add logcounts for viz
message("normalizing counts - ", Sys.time())
set.seed(1000)
sce <- computeSumFactors(sce, cluster=sce$k_50_label)
sce <- logNormCounts(sce)

message("saving data - ", Sys.time())
save(sce, here("snRNAseq_hpc","processed-data", "sce", "sce_clustered.rda"))

pdf(here("snRNAseq_hpc","plots","markers_k50.pdf"))
plotExpression(sce,features=c('SYT1','SLC17A7','SLC17A6',
                              'GAD1','GAD2','MBP',
                              'GFAP','TTR','CSF1R',
                              'PDGFRA','FLT1','TNNT2'),
               x="k_50_label", colour_by="k_50_label", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", 
               width = 0.3)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()