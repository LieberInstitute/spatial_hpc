library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("here")
library("sessioninfo")
library("batchelor")
library("BiocSingular")
library('ggrastr')

##load sce
load(here::here("snRNAseq_hpc","processed-data", "sce", "sce_post_qc.rda"))
dim(sce)
##load palettes
load(here::here('plots','snRNAseq_palettes.rda'))
# [1] 36601 86905

##make sure these are factors
sce$brnum<-factor(sce$brnum)
sce$round<-factor(sce$round)
sce$sort<-factor(sce$sort)

##make sure feature names unique
rownames(sce)<-uniquifyFeatureNames(rowData(sce)$gene_id,rowData(sce)$gene_name)

#Remove nuclei from problematic sample
sce$discard<-ifelse(sce$sample_ID %in% "17c-scp" & sce$detected < 5000,T,F)
sce<-sce[,sce$discard==F]
dim(sce)
# [1] 36601 80594

##feature selection using deviance w/ poisson model. Correcting for batch
##We'll use brnum as the batch variable here
set.seed(800)
sce <- devianceFeatureSelection(sce,
                                assay = "counts", fam = "poisson",
                                sorted = T,batch=sce$brnum)

pdf(here("snRNAseq_hpc","plots", "build_sce", "binomial_deviance.pdf"))
plot(sort(rowData(sce)$poisson_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylim=c(0,1000000),
     ylab = "poisson deviance", main = "Feature Selection with Deviance"
)

abline(v = 2000, lty = 2, col = "blue")
dev.off()

hdg<-rownames(counts(sce))[1:2000]
res<-sce[rownames(counts(sce)) %in% hdg,]
set.seed(913)

res <- nullResiduals(res,
                      fam = "poisson",
                          type = "pearson")


set.seed(915)
message("running PCA - ", Sys.time())
res <- scater::runPCA(res,ntop=2000,
                      exprs_values='poisson_pearson_residuals',
                      ncomponents=100,
                      BSPARAM=IrlbaParam())

reducedDim(sce,'PCA')<-reducedDim(res,'PCA')

set.seed(125)
message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "PCA")

# Run mnn
message("running MNN - ", Sys.time())
set.seed(1788)
mnn<-batchelor::reducedMNN(reducedDim(sce,'PCA'),batch=sce$brnum,k=50)
reducedDim(sce,'MNN')<-mnn$corrected
set.seed(1788)
mnn<-batchelor::reducedMNN(reducedDim(sce,'MNN'),batch=sce$round,k=50)
reducedDim(sce,'MNN')<-mnn$corrected

set.seed(2351)
message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "MNN")



message("normalizing counts - ", Sys.time())
set.seed(100)
clust <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clust, min.mean=0.1)
sce <- logNormCounts(sce)

##visualize neuronal marker gene
plotUMAP(sce,colour_by='SYT1',point_size=0.25)

##Initial clustering
message("running buildSNNGraph - ", Sys.time())
snn.gr <- buildSNNGraph(sce, k = 5, use.dimred = "MNN",type='jaccard')

set.seed(100)
clust_5 <- igraph::cluster_louvain(snn.gr,resolution=1.25)$membership
table(clust_5)
sce$k_5_louvain_initial<-factor(clust_5)


##flag low-quality clusters (see QC_snRNAseq_suppFigure_plots.R and supp fig 9)
sce$discard<-ifelse(sce$k_5_louvain_initial %in% c(6,8,37,46),T,F)
sce<-sce[,sce$discard==F]

save(sce,file=here::here("snRNAseq_hpc",'processed-data','sce','sce_clustered_round2.rda'))
##########Rerun all steps###############
##feature selection using deviance w/ poisson model. Correcting for batch
##We'll use brnum as the batch variable for feature selection and dimred
##Large differences in anatomy
set.seed(800)
sce <- devianceFeatureSelection(sce,
                                assay = "counts", fam = "poisson",
                                sorted = T,batch=sce$brnum)

hdg<-rownames(counts(sce))[1:2000]
res<-sce[rownames(counts(sce)) %in% hdg,]
set.seed(913)
res <- nullResiduals(res,
                     fam = "poisson",
                     type = "pearson")

set.seed(915)
message("running PCA - ", Sys.time())
res <- scater::runPCA(res,ntop=2000,
                      exprs_values='poisson_pearson_residuals',
                      ncomponents=100,
                      BSPARAM=IrlbaParam())

reducedDim(sce,'PCA')<-reducedDim(res,'PCA')

set.seed(125)
message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "PCA")

# Run mnn
message("running MNN - ", Sys.time())
set.seed(1788)
mnn<-batchelor::reducedMNN(reducedDim(sce,'PCA'),batch=sce$brnum,k=50)
reducedDim(sce,'MNN')<-mnn$corrected
set.seed(1788)
mnn<-batchelor::reducedMNN(reducedDim(sce,'MNN'),batch=sce$round,k=50)
reducedDim(sce,'MNN')<-mnn$corrected

set.seed(2351)
message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "MNN")

##Initial clustering
message("running buildSNNGraph - ", Sys.time())
snn.gr <- buildSNNGraph(sce, k = 5, use.dimred = "MNN",type='jaccard')

set.seed(100)
clust_5 <- igraph::cluster_louvain(snn.gr,resolution=1.5)$membership
table(clust_5)
sce$k_5_louvain2<-factor(clust_5)

#Flag
sce$discard2<-ifelse(sce$k_5_louvain %in% c(17,23),T,F)
save(sce,file=here::here("snRNAseq_hpc",'processed-data','sce',
                         'sce_clustered_round3.rda'))
sce<-sce[,sce$discard2==F]

set.seed(2351)
message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "MNN")

save(sce,file=here::here("snRNAseq_hpc",'processed-data','sce',
                         'sce_clustered_round4.rda'))

tab<-data.frame('cluster'=1:62,'mid'=rep(NA,62),
                'fine'=rep(NA,62),
                'type'=rep(NA,62))
########Neurons
##Granule cells
tab$mid[tab$cluster %in% c(33,12,55,31,49)]<-"ExcN"
tab$fine[tab$cluster %in% c(33,12,55,31,49)]<-'GC'
tab$type[tab$cluster %in% c(33,12,55,31,49)]<-paste0('GC.',c(1:5))

##Mossy cells
tab$mid[tab$cluster %in% c(21)]<-"ExcN"
tab$fine[tab$cluster %in% c(21)]<-'CA2-4'
tab$type[tab$cluster %in% c(21)]<-'MC'
##CA3
tab$mid[tab$cluster %in% c(11,38)]<-"ExcN"
tab$fine[tab$cluster %in% c(11,38)]<-'CA2-4'
tab$type[tab$cluster %in% c(11,38)]<-paste0('CA3.',c(1:2))
##CA2
tab$mid[tab$cluster %in% c(45)]<-"ExcN"
tab$fine[tab$cluster %in% c(45)]<-'CA2-4'
tab$type[tab$cluster %in% c(45)]<-'CA2'

##CA1
tab$mid[tab$cluster %in% c(26)]<-"ExcN"
tab$fine[tab$cluster %in% c(26)]<-'CA1/ProS'
tab$type[tab$cluster %in% c(26)]<-'CA1'

##ProS
tab$mid[tab$cluster %in% c(14)]<-"ExcN"
tab$fine[tab$cluster %in% c(14)]<-'CA1/ProS'
tab$type[tab$cluster %in% c(14)]<-'ProS'

##Sub_sup
tab$mid[tab$cluster %in% c(25)]<-"ExcN"
tab$fine[tab$cluster %in% c(25)]<-'Sub.1'
tab$type[tab$cluster %in% c(25)]<-'Sub.1'
##Sub_deep
tab$mid[tab$cluster %in% c(22)]<-"ExcN"
tab$fine[tab$cluster %in% c(22)]<-'Sub.2'
tab$type[tab$cluster %in% c(22)]<-'Sub.2'

##RHP
##L2/3
tab$mid[tab$cluster %in% c(13,51,54,28,37,50)]<-"ExcN"
tab$fine[tab$cluster %in% c(13,51,54,28,37,50)]<-'L2/3'
tab$type[tab$cluster %in% c(13,51,54,28,37,50)]<-paste0('L2/3.',c(1:6))

##L5
tab$mid[tab$cluster %in% c(5,43)]<-"ExcN"
tab$fine[tab$cluster %in% c(5,43)]<-'L5/6'
tab$type[tab$cluster %in% c(5,43)]<-paste0('L5/6.',c(1,2))

##L6
tab$mid[tab$cluster %in% c(32,15)]<-"ExcN"
tab$fine[tab$cluster %in% c(32,15)]<-'L6/6b'
tab$type[tab$cluster %in% c(32,15)]<-paste0('L6.',c(1,2))

##L6b
tab$mid[tab$cluster %in% c(35)]<-"ExcN"
tab$fine[tab$cluster %in% c(35)]<-'L6/6b'
tab$type[tab$cluster %in% c(35)]<-'L6b'

##Thalamus
tab$mid[tab$cluster %in% c(4)]<-"ExcN"
tab$fine[tab$cluster %in% c(4)]<-'Thal'
tab$type[tab$cluster %in% c(4)]<-'Thal'



##AHi
tab$mid[tab$cluster %in% c(57,60,46,61)]<-"ExcN"
tab$fine[tab$cluster %in% c(57,60,46,61)]<-'Amy'
tab$type[tab$cluster %in% c(57,60,46,61)]<-paste0('AHi.',c(1:4))
##HATA
tab$mid[tab$cluster %in% c(53)]<-"ExcN"
tab$fine[tab$cluster %in% c(53)]<-'HATA'
tab$type[tab$cluster %in% c(53)]<-'HATA'

##Cajal Retzius
tab$mid[tab$cluster %in% c(58)]<-"ExcN"
tab$fine[tab$cluster %in% c(58)]<-'Cajal'
tab$type[tab$cluster %in% c(58)]<-'Cajal'

##GABAergics
tab$mid[tab$cluster %in% c(59)]<-"InhN"
tab$fine[tab$cluster %in% c(59)]<-'GABA.PENK'
tab$type[tab$cluster %in% c(59)]<-'PENK'

tab$mid[tab$cluster %in% c(10)]<-"InhN"
tab$fine[tab$cluster %in% c(10)]<-'GABA.MGE'
tab$type[tab$cluster %in% c(10)]<-'SST'

tab$mid[tab$cluster %in% c(52)]<-"InhN"
tab$fine[tab$cluster %in% c(52)]<-'GABA.MGE'
tab$type[tab$cluster %in% c(52)]<-'CORT'

tab$mid[tab$cluster %in% c(16)]<-"InhN"
tab$fine[tab$cluster %in% c(16)]<-'GABA.MGE'
tab$type[tab$cluster %in% c(16)]<-'PV.FS'

tab$mid[tab$cluster %in% c(56)]<-"InhN"
tab$fine[tab$cluster %in% c(56)]<-'GABA.MGE'
tab$type[tab$cluster %in% c(56)]<-'CRABP1'

tab$mid[tab$cluster %in% c(7)]<-"InhN"
tab$fine[tab$cluster %in% c(7)]<-'GABA.MGE'
tab$type[tab$cluster %in% c(7)]<-'C1QL1'

tab$mid[tab$cluster %in% c(3)]<-"InhN"
tab$fine[tab$cluster %in% c(3)]<-'GABA.LAMP5'
tab$type[tab$cluster %in% c(3)]<-'LAMP5.MGE'

tab$mid[tab$cluster %in% c(18)]<-"InhN"
tab$fine[tab$cluster %in% c(18)]<-'GABA.LAMP5'
tab$type[tab$cluster %in% c(18)]<-'LAMP5.CGE'

tab$mid[tab$cluster %in% c(39)]<-"InhN"
tab$fine[tab$cluster %in% c(39)]<-'GABA.LAMP5'
tab$type[tab$cluster %in% c(39)]<-'CXCL14'

tab$mid[tab$cluster %in% c(40)]<-"InhN"
tab$fine[tab$cluster %in% c(40)]<-'GABA.CGE'
tab$type[tab$cluster %in% c(40)]<-'HTR3A'

tab$mid[tab$cluster %in% c(29)]<-"InhN"
tab$fine[tab$cluster %in% c(29)]<-'GABA.CGE'
tab$type[tab$cluster %in% c(29)]<-'VIP'

######NNC
tab$mid[tab$cluster %in% c(1,9)]<-'Oligo'
tab$fine[tab$cluster %in% c(1,9)]<-'Oligo'
tab$type[tab$cluster %in% c(1,9)]<-paste0('Oligo.',c(1:2))

tab$mid[tab$cluster %in% c(8,20,27)]<-'Astro'
tab$fine[tab$cluster %in% c(8,20,27)]<-'Astro'
tab$type[tab$cluster %in% c(8,20,27)]<-paste0('Astro.',c(1:3))

tab$mid[tab$cluster %in% c(24)]<-'OPC'
tab$fine[tab$cluster %in% c(24)]<-'OPC'
tab$type[tab$cluster %in% c(24)]<-'OPC'

tab$mid[tab$cluster %in% c(44)]<-'OPC'
tab$fine[tab$cluster %in% c(44)]<-'OPC'
tab$type[tab$cluster %in% c(44)]<-'COP'

tab$mid[tab$cluster %in% c(2,6)]<-'Micro/Macro/T'
tab$fine[tab$cluster %in% c(2,6)]<-'Micro/Macro/T'
tab$type[tab$cluster %in% c(2,6)]<-paste0('Micro.',c(1:2))

tab$mid[tab$cluster %in% c(42)]<-'Micro/Macro/T'
tab$fine[tab$cluster %in% c(42)]<-'Micro/Macro/T'
tab$type[tab$cluster %in% c(42)]<-'Macro/Tcell'

tab$mid[tab$cluster %in% c(48)]<-'CSF'
tab$fine[tab$cluster %in% c(48)]<-'Ependy'
tab$type[tab$cluster %in% c(48)]<-'Ependy'

tab$mid[tab$cluster %in% c(34,36,62)]<-'CSF'
tab$fine[tab$cluster %in% c(34,36,62)]<-'Choroid'
tab$type[tab$cluster %in% c(34,36,62)]<-paste0('CP.',c(1:3))

tab$mid[tab$cluster %in% c(19,47)]<-'Vascular'
tab$fine[tab$cluster %in% c(19,47)]<-'Vascular'
tab$type[tab$cluster %in% c(19,47)]<-paste0('Endo.',c(1:2))

tab$mid[tab$cluster %in% c(30)]<-'Vascular'
tab$fine[tab$cluster %in% c(30)]<-'Vascular'
tab$type[tab$cluster %in% c(30)]<-'PC/SMC'

tab$mid[tab$cluster %in% c(41)]<-'Vascular'
tab$fine[tab$cluster %in% c(41)]<-'Vascular'
tab$type[tab$cluster %in% c(41)]<-'VLMC'

##split up the L2/3s


sce$mid.cell.class<-factor(tab$mid[match(sce$k_5_louvain,tab$cluster)],
                           levels=c('ExcN','InhN','Micro/Macro/T',
                                    'Astro','Oligo','OPC','CSF','Vascular'))
sce$fine.cell.class<-factor(tab$fine[match(sce$k_5_louvain,tab$cluster)],
                            levels=c('GC','CA2-4','CA1/ProS','Sub.1','Sub.2',
                                     'L6/6b','L5/6','L2/3',
                                     'HATA','Amy','Thal','Cajal',
                                     'GABA.PENK','GABA.MGE','GABA.LAMP5','GABA.CGE',
                                     'Micro/Macro/T','Astro',
                                     'Oligo','OPC','Ependy','Choroid','Vascular'))
sce$superfine.cell.class<-factor(tab$type[match(sce$k_5_louvain,tab$cluster)],
                      levels=c("GC.1", "GC.2", "GC.3", "GC.4", "GC.5", "MC",
                               "CA3.1", "CA3.2", "CA2", "CA1", "ProS", "Sub.1",
                               "Sub.2",
                               "L6.2", "L6.1", "L6b",
                               "L5.1", "L5.2",
                               "L2/3.1","L2/3.5", "L2/3.2", "L2/3.3", "L2/3.4",
                               "L2/3.6",
                               "HATA", "AHi.1", "AHi.2","AHi.3","AHi.4","Thal",
                               "Cajal", "PENK", "SST", "CORT", "PV.FS",
                               "CRABP1", "C1QL1", "LAMP5.MGE", "LAMP5.CGE",
                               "CXCL14", "HTR3A", "VIP",
                               "Micro.1", "Micro.2", "Macro/Tcell",
                               "Astro.1", "Astro.2", "Astro.3",
                               "Oligo.1", "Oligo.2",
                               "OPC", "COP",
                               "Ependy",
                               "CP.1", "CP.2", "CP.3",
                               "Endo.1", "Endo.2",
                               "PC/SMC", "VLMC"))

sce$broad.cell.class<-factor(
    ifelse(sce$mid.cell.class %in% levels(sce$mid.cell.class)[1:2],'Neuron',
           ifelse(sce$mid.cell.class %in% levels(sce$mid.cell.class)[3],'Micro',
                  ifelse(sce$mid.cell.class %in% levels(sce$mid.cell.class)[4],'Astro',
                         ifelse(sce$mid.cell.class %in% levels(sce$mid.cell.class)[5],'Oligo',
                                "Other")))),
    levels=c('Neuron','Micro',
             'Astro','Oligo','Other')
)

###split up the layer 2/3s
sce$fine.cell.class<-factor(ifelse(sce$superfine.cell.class %in% c('L2/3.1','L2/3.5'),'L2/3.PrS.PaS',
                            ifelse(sce$fine.type %in% c('L2/3.2','L2/3.4','L2/3.6','L2/3.3'),'L2/3.PrS.Ent',
                            as.character(sce$fine.cellclass))),
                            levels=c('GC','CA2-4','CA1/ProS','Sub.1','Sub.2',
                                     'L6/6b','L5/6','L2/3.PrS.PaS','L2/3.PrS.Ent',
                                     'HATA','Amy','Thal','Cajal',
                                     'GABA.PENK','GABA.MGE','GABA.LAMP5','GABA.CGE',
                                     'Micro/Macro/T','Astro',
                                     'Oligo','OPC','Ependy','Choroid','Vascular'))

save(sce,file=here::here("snRNAseq_hpc",'processed-data','sce',
                         'sce_final.rda'))
