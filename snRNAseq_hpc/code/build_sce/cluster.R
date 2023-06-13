library("SingleCellExperiment")
library("jaffelab")
library("scater")
library("scran")
library("here")
library("sessioninfo")

load(here("snRNAseq_hpc","processed-data", "sce", "sce_MNN.rda"), verbose = TRUE)

message("running buildSNNGraph - ", Sys.time())
snn.gr <- buildSNNGraph(sce, k = 10, use.dimred = "MNN",type='jaccard')


message("running louvain - ", Sys.time())
set.seed(100)
clust_5 <- igraph::cluster_louvain(snn.gr,resolution=1.5)$membership
table(clust_5)
sce$k_10_louvain<-factor(clust_5)

sce$discard<-ifelse(sce$k_10_louvain %in% c(6,7,48),T,F)
sce<-sce[,sce$discard==F]

set.seed(000)
clust_5 <- igraph::cluster_louvain(snn.gr,resolution=2)$membership
table(clust_5)
sce$k_5_louvain_4<-factor(clust_5)


table(clust_50)
table(clust_5)

##add to sce
sce$k_5_louvain_4<-factor(clust_5)
sce$k_25_pca_walk<-clust_25

##make some prelim plots
pdf(here("snRNAseq_hpc","plots","UMAP_k20.pdf"))
plotUMAP(sce,colour_by='k_20_label',text_by='k_20_label')
dev.off()

pdf(here("snRNAseq_hpc","plots","UMAP_k50.pdf"))
plotUMAP(sce,colour_by='k_50_label',text_by='k_50_label')
dev.off()

message("saving data - ", Sys.time())
save(sce, file=here("snRNAseq_hpc","processed-data", "sce", "sce_clustered.rda"))

##add logcounts for viz
message("normalizing counts - ", Sys.time())
set.seed(1000)
sce <- computeSumFactors(sce, cluster=sce$k_50_label)
sce <- logNormCounts(sce)

message("saving data - ", Sys.time())
save(sce, file=here("snRNAseq_hpc","processed-data", "sce", "sce_clustered.rda"))

pdf(here("snRNAseq_hpc","plots","markers_k50.pdf"))
plotExpression(sce,features=c('SYT1','SLC17A7','SLC17A6',
                              'GAD1','GAD2','MBP',
                              'GFAP','TTR','CSF1R',
                              'PDGFRA','FLT1','TNNT2'),
               x="k_50_label", colour_by="k_50_label", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
               width = 0.3)
dev.off()

x<-RcppML::nmf(
    assay(sce,'logcounts'),
    75,
    tol = 1e-04,
    maxit = 100,
    verbose = TRUE,
    L1 = c(0, 0),
    mask_zeros = FALSE,
    diag = TRUE,
    nonneg = TRUE,
    seed=122
)

####First remove low quality neuron clusters from first run and rerun all feature selection/dimred/batch correction steps
message("running louvain - ", Sys.time())
set.seed(100)
clust_5 <- igraph::cluster_louvain(snn.gr,resolution=1.2)$membership
table(clust_5)
sce$cluster_harmony<-factor(clust_5)



tab<-data.frame('cluster'=1:51,'annotation'=rep(NA,51),'annoType'=rep(NA,51))
########Neurons
##Granule cells
tab$annotation[tab$cluster %in% c(10,28,39,40)]<-'EXC'
tab$annoType[tab$cluster %in% c(10,28,39,40)]<-paste0('GC.',c(1:4))
##Mossy cells
tab$annotation[tab$cluster %in% c(18)]<-'EXC'
tab$annoType[tab$cluster %in% c(18)]<-'MC'
##CA3
tab$annotation[tab$cluster %in% c(9)]<-'EXC'
tab$annoType[tab$cluster %in% c(9)]<-'CA3'
##CA3
tab$annotation[tab$cluster %in% c(33)]<-'EXC'
tab$annoType[tab$cluster %in% c(33)]<-'CA3.deep'
##CA2
tab$annotation[tab$cluster %in% c(36)]<-'EXC'
tab$annoType[tab$cluster %in% c(36)]<-'CA2'
##CA1
tab$annotation[tab$cluster %in% c(11)]<-'EXC'
tab$annoType[tab$cluster %in% c(11)]<-'CA1'

##Sub_sup
tab$annotation[tab$cluster %in% c(23)]<-'EXC'
tab$annoType[tab$cluster %in% c(23)]<-'Sub.Sup'
##Sub_deep
tab$annotation[tab$cluster %in% c(19)]<-'EXC'
tab$annoType[tab$cluster %in% c(19)]<-'Sub.Deep'

##RHP
tab$annotation[tab$cluster %in% c(25,29,38,43,45)]<-'EXC'
tab$annoType[tab$cluster %in% c(25,29,38,43,45)]<-paste0('RHP.L2.',c(1:5))

tab$annotation[tab$cluster %in% c(5)]<-'EXC'
tab$annoType[tab$cluster %in% c(5)]<-'RHP.L3'

tab$annotation[tab$cluster %in% c(25,29,38,43,45)]<-'EXC'
tab$annoType[tab$cluster %in% c(25,29,38,43,45)]<-paste0('RHP.L2.',c(1:5))

tab$annotation[tab$cluster %in% c(21)]<-'EXC'
tab$annoType[tab$cluster %in% c(21)]<-'RHP.L5'

tab$annotation[tab$cluster %in% c(12,14)]<-'EXC'
tab$annoType[tab$cluster %in% c(12,14)]<-paste0('RHP.L6.',c(1:2))

tab$annotation[tab$cluster %in% c(20)]<-'EXC'
tab$annoType[tab$cluster %in% c(20)]<-paste0('RHP.L6b')

tab$annotation[tab$cluster %in% c(4)]<-'EXC'
tab$annoType[tab$cluster %in% c(4)]<-paste0('Thal')



##Amygdala
tab$annotation[tab$cluster %in% c(41,42,48)]<-'EXC'
tab$annoType[tab$cluster %in% c(41,42,48)]<-paste0('Amy.',c(1:3))
##HATA
tab$annotation[tab$cluster %in% c(32)]<-'EXC'
tab$annoType[tab$cluster %in% c(32)]<-'HATA'

##Cajal Retzius
tab$annotation[tab$cluster %in% c(46)]<-'EXC'
tab$annoType[tab$cluster %in% c(46)]<-'Cajal'

##GABAergics
tab$annotation[tab$cluster %in% c(47)]<-'INH'
tab$annoType[tab$cluster %in% c(47)]<-'MEIS2'

tab$annotation[tab$cluster %in% c(13)]<-'INH'
tab$annoType[tab$cluster %in% c(13)]<-'PVALB'

tab$annotation[tab$cluster %in% c(8)]<-'INH'
tab$annoType[tab$cluster %in% c(8)]<-'SST'

tab$annotation[tab$cluster %in% c(49)]<-'INH'
tab$annoType[tab$cluster %in% c(49)]<-'NKX2-1'

tab$annotation[tab$cluster %in% c(3)]<-'INH'
tab$annoType[tab$cluster %in% c(3)]<-'LAMP5 LHX6'

tab$annotation[tab$cluster %in% c(15)]<-'INH'
tab$annoType[tab$cluster %in% c(15)]<-'LAMP5 KIT'

tab$annotation[tab$cluster %in% c(34)]<-'INH'
tab$annoType[tab$cluster %in% c(34)]<-'HTR3A'

tab$annotation[tab$cluster %in% c(26)]<-'INH'
tab$annoType[tab$cluster %in% c(26)]<-'VIP'

tab$annotation[tab$cluster %in% c(44)]<-'INH'
tab$annoType[tab$cluster %in% c(44)]<-'CORT'

######NNC
tab$annotation[tab$cluster %in% c(51,30)]<-'NNC'
tab$annoType[tab$cluster %in% c(51,30)]<-paste0('CP.',c(1:2))

tab$annotation[tab$cluster %in% c(1,7)]<-'NNC'
tab$annoType[tab$cluster %in% c(1,7)]<-paste0('Oligo.',c(1:2))

tab$annotation[tab$cluster %in% c(31)]<-'NNC'
tab$annoType[tab$cluster %in% c(31)]<-'Oligo.NB'

tab$annotation[tab$cluster %in% c(50)]<-'NNC'
tab$annoType[tab$cluster %in% c(50)]<-'COP'

tab$annotation[tab$cluster %in% c(22)]<-'NNC'
tab$annoType[tab$cluster %in% c(22)]<-'OPC'

tab$annotation[tab$cluster %in% c(6,17)]<-'NNC'
tab$annoType[tab$cluster %in% c(6,17)]<-paste0('Astro.',c(1:2))

tab$annotation[tab$cluster %in% c(37)]<-'NNC'
tab$annoType[tab$cluster %in% c(37)]<-'Ependy'

tab$annotation[tab$cluster %in% c(37)]<-'NNC'
tab$annoType[tab$cluster %in% c(37)]<-'Ependy'

tab$annotation[tab$cluster %in% c(2)]<-'NNC'
tab$annoType[tab$cluster %in% c(2)]<-'Micro'

tab$annotation[tab$cluster %in% c(24)]<-'NNC'
tab$annoType[tab$cluster %in% c(24)]<-'Tcell/Macro'

tab$annotation[tab$cluster %in% c(27)]<-'NNC'
tab$annoType[tab$cluster %in% c(27)]<-'SMC'

tab$annotation[tab$cluster %in% c(16)]<-'NNC'
tab$annoType[tab$cluster %in% c(16)]<-'Endo/PC'

tab$annotation[tab$cluster %in% c(35)]<-'NNC'
tab$annoType[tab$cluster %in% c(35)]<-'VLMC'

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
