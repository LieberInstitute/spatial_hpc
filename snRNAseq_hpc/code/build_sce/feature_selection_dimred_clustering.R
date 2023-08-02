library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("here")
library("sessioninfo")
library("harmony")
library("BiocSingular")
library('ggrastr')

## Load sce
load(here("snRNAseq_hpc","processed-data", "sce", "sce_post_qc.rda"))
dim(sce)
#
# [1] 36601 86905
sce$brnum<-factor(tab$Brain[match(sce$Sample,tab[,1])])
sce$round<-factor(tab$Round[match(sce$Sample,tab[,1])])
sce$sort<-factor(tab$PI.NeuN[match(sce$Sample,tab[,1])])
#Remove nuclei from problematic sample
sce$discard<-ifelse(sce$Sample %in% "17c-scp" & sce$detected < 5000,T,F)
sce<-sce[,sce$discard==F]
# [1] 36601 80594
##feature selection using deviance w/ poisson model. Correcting for batch
##We'll use brnum as the batch variable for feature selection and dimred
##Large differences in anatomy
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
# # Initialize an empty list to store the results.
# res_list <- list()
#
# # Loop over the subsets defined by brnum.
# for (i in seq_along(splitit(res$brnum))) {
#     # Subset the sce object.
#     res_temp <- res[, splitit(res$brnum)[[i]]]
#
#     message("running nullResiduals - ", Sys.time())
#
#     # Apply the nullResiduals function.
#     res_temp <- nullResiduals(res_temp,
#                               fam = "poisson",
#                               type = "pearson"
#     )
#
#     # Add the result to the list.
#     res_list[[i]] <- res_temp
# }
#
# # Combine the results back into a single sce object.
# res_combined <- do.call(cbind, res_list)
# idx <- match(colnames(res), colnames(res_combined))
#
# # Use the index to order the columns of 'res_combined'
# res_combined <- res_combined[,idx]
#
# # Check that the column names are the same as in the original res object.
# all(colnames(res) == colnames(res_combined))
#
# res<-res_combined
# rm(res_combined)

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

pdf(here("snRNAseq_hpc","plots","UMAPs_unccorrected.pdf"))
plotUMAP(sce,colour_by='Sample',point_size=0.25)
plotUMAP(sce,colour_by='brnum',point_size=0.25)
plotUMAP(sce,colour_by='sort',point_size=0.25)
plotUMAP(sce,colour_by='round',point_size=0.25)
dev.off()

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
plotUMAP(sce,colour_by='SYT1',point_size=0.25)


message("normalizing counts - ", Sys.time())
set.seed(100)
clust <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clust, min.mean=0.1)
sce <- logNormCounts(sce)

pdf(here("snRNAseq_hpc","plots","UMAPs_corrected.pdf"))
plotUMAP(sce,colour_by='Sample',point_size=0.25)
plotUMAP(sce,colour_by='brnum',point_size=0.25)
plotUMAP(sce,colour_by='sort',point_size=0.25)
plotUMAP(sce,colour_by='round',point_size=0.25)
plotUMAP(sce,colour_by='SYT1',point_size=0.25)
plotUMAP(sce,colour_by='detected',point_size=0.25)
dev.off()


##Initial clustering
message("running buildSNNGraph - ", Sys.time())
snn.gr <- buildSNNGraph(sce, k = 5, use.dimred = "MNN",type='jaccard')

set.seed(100)
clust_5 <- igraph::cluster_louvain(snn.gr,resolution=1.25)$membership
table(clust_5)
sce$k_5_louvain_initial<-factor(clust_5)

#pdf(here("snRNAseq_hpc","plots","UMAPs_corrected.pdf"))
plotUMAP(sce,text_by='k_5_louvain_initial',
         colour_by='k_5_louvain_initial',point_size=0.25)



##flag some low-quality neuron clusters
sce$discard<-ifelse(sce$k_5_louvain_initial %in% c(6,8,18,25,39),T,F)
sce<-sce[,sce$discard==F]

##make neuron vs non-neuron comparisons
#index<-c(2,9,21,25,10,1,40,36,35,39,53,27,20,29)
#sce$neuron<-ifelse(sce$k_5_louvain_initial3 %in% index,F,T)


##Rerun all steps
##feature selection using deviance w/ poisson model. Correcting for batch
##We'll use brnum as the batch variable for feature selection and dimred
##Large differences in anatomy
set.seed(800)
sce <- devianceFeatureSelection(sce,
                                assay = "counts", fam = "poisson",
                                sorted = T,batch=sce$brnum)

pdf(here("snRNAseq_hpc","plots", "build_sce", "binomial_deviance_corrected.pdf"))
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

pdf(here("snRNAseq_hpc","plots","UMAPs_unccorrected_filtered.pdf"))
plotUMAP(sce,colour_by='Sample',point_size=0.25)
plotUMAP(sce,colour_by='brnum',point_size=0.25)
plotUMAP(sce,colour_by='sort',point_size=0.25)
plotUMAP(sce,colour_by='round',point_size=0.25)
dev.off()

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

pdf(here("snRNAseq_hpc","plots","UMAPs_corrected.pdf"))
plotUMAP(sce,colour_by='Sample',point_size=0.25)
plotUMAP(sce,colour_by='brnum',point_size=0.25)
plotUMAP(sce,colour_by='sort',point_size=0.25)
plotUMAP(sce,colour_by='round',point_size=0.25)
plotUMAP(sce,colour_by='SYT1',point_size=0.25)
plotUMAP(sce,colour_by='detected',point_size=0.25)
dev.off()


##Initial clustering
message("running buildSNNGraph - ", Sys.time())
snn.gr <- buildSNNGraph(sce, k = 5, use.dimred = "MNN",type='jaccard')

set.seed(100)
clust_5 <- igraph::cluster_louvain(snn.gr,resolution=1.5)$membership
table(clust_5)
sce$k_5_louvain2<-factor(clust_5)

#Flag
sce$discard2<-ifelse(sce$k_5_louvain %in% c(17,23),T,F)
sce<-sce[,sce$discard2==F]

set.seed(2351)
message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "MNN")
marks<-findMarkers(sce,groups=sce$k_5_louvain,pval.type='all',direction='up',test.type='wilcox')

tab<-data.frame('cluster'=1:62,'broad'=rep(NA,62),
                'anno'=rep(NA,62),
                'type'=rep(NA,62))
########Neurons
##Granule cells
tab$broad[tab$cluster %in% c(33,12,55,31,49)]<-'EXC'
tab$anno[tab$cluster %in% c(33,12,55,31,49)]<-'GC'
tab$type[tab$cluster %in% c(33,12,55,31,49)]<-paste0('GC.',c(1:5))

##Mossy cells
tab$broad[tab$cluster %in% c(21)]<-'EXC'
tab$anno[tab$cluster %in% c(21)]<-'CA2-4'
tab$type[tab$cluster %in% c(21)]<-'MC'
##CA3
tab$broad[tab$cluster %in% c(11,38)]<-'EXC'
tab$anno[tab$cluster %in% c(11,38)]<-'CA2-4'
tab$type[tab$cluster %in% c(11,38)]<-paste0('CA3.',c(1:2))
##CA2
tab$broad[tab$cluster %in% c(45)]<-'EXC'
tab$anno[tab$cluster %in% c(45)]<-'CA2-4'
tab$type[tab$cluster %in% c(45)]<-'CA2'

##CA1
tab$broad[tab$cluster %in% c(26)]<-'EXC'
tab$anno[tab$cluster %in% c(26)]<-'CA1'
tab$type[tab$cluster %in% c(26)]<-'CA1'

##ProS
tab$broad[tab$cluster %in% c(14)]<-'EXC'
tab$anno[tab$cluster %in% c(14)]<-'ProS/Sub'
tab$type[tab$cluster %in% c(14)]<-'ProS'

##Sub_sup
tab$broad[tab$cluster %in% c(25)]<-'EXC'
tab$anno[tab$cluster %in% c(25)]<-'ProS/Sub'
tab$type[tab$cluster %in% c(25)]<-'Sub.1'
##Sub_deep
tab$broad[tab$cluster %in% c(22)]<-'EXC'
tab$anno[tab$cluster %in% c(22)]<-'ProS/Sub'
tab$type[tab$cluster %in% c(22)]<-'Sub.2'

##RHP
##L2/3
tab$broad[tab$cluster %in% c(13,51,54,28,37,50)]<-'EXC'
tab$anno[tab$cluster %in% c(13,51,54,28,37,50)]<-'L2/3'
tab$type[tab$cluster %in% c(13,51,54,28,37,50)]<-paste0('L2/3.',c(1:6))

##L5
tab$broad[tab$cluster %in% c(5,43)]<-'EXC'
tab$anno[tab$cluster %in% c(5,43)]<-'L5'
tab$type[tab$cluster %in% c(5,43)]<-paste0('L5.',c(1,2))

##L6
tab$broad[tab$cluster %in% c(32,15)]<-'EXC'
tab$anno[tab$cluster %in% c(32,15)]<-'L6/6b'
tab$type[tab$cluster %in% c(32,15)]<-paste0('L6.',c(1,2))

##L6b
tab$broad[tab$cluster %in% c(35)]<-'EXC'
tab$anno[tab$cluster %in% c(35)]<-'L6/6b'
tab$type[tab$cluster %in% c(35)]<-'L6b'

##Thalamus
tab$broad[tab$cluster %in% c(4)]<-'EXC'
tab$anno[tab$cluster %in% c(4)]<-'Thal'
tab$type[tab$cluster %in% c(4)]<-'Thal'



##AHi
tab$broad[tab$cluster %in% c(57,60,46,61)]<-'EXC'
tab$anno[tab$cluster %in% c(57,60,46,61)]<-'HATA/AHi'
tab$type[tab$cluster %in% c(57,60,46,61)]<-paste0('AHi.',c(1:4))
##HATA
tab$broad[tab$cluster %in% c(53)]<-'EXC'
tab$anno[tab$cluster %in% c(53)]<-'HATA/AHi'
tab$type[tab$cluster %in% c(53)]<-'HATA'

##Cajal Retzius
tab$broad[tab$cluster %in% c(58)]<-'EXC'
tab$anno[tab$cluster %in% c(58)]<-'Cajal'
tab$type[tab$cluster %in% c(58)]<-'Cajal'

##GABAergics
tab$broad[tab$cluster %in% c(59)]<-'INH'
tab$anno[tab$cluster %in% c(59)]<-'GABA'
tab$type[tab$cluster %in% c(59)]<-'PENK'

tab$broad[tab$cluster %in% c(10)]<-'INH'
tab$anno[tab$cluster %in% c(10)]<-'GABA'
tab$type[tab$cluster %in% c(10)]<-'SST'

tab$broad[tab$cluster %in% c(52)]<-'INH'
tab$anno[tab$cluster %in% c(52)]<-'GABA'
tab$type[tab$cluster %in% c(52)]<-'CORT'

tab$broad[tab$cluster %in% c(16)]<-'INH'
tab$anno[tab$cluster %in% c(16)]<-'GABA'
tab$type[tab$cluster %in% c(16)]<-'PV.FS'

tab$broad[tab$cluster %in% c(56)]<-'INH'
tab$anno[tab$cluster %in% c(56)]<-'GABA'
tab$type[tab$cluster %in% c(56)]<-'CRABP1'

tab$broad[tab$cluster %in% c(7)]<-'INH'
tab$anno[tab$cluster %in% c(7)]<-'GABA'
tab$type[tab$cluster %in% c(7)]<-'C1QL1'

tab$broad[tab$cluster %in% c(3)]<-'INH'
tab$anno[tab$cluster %in% c(3)]<-'GABA'
tab$type[tab$cluster %in% c(3)]<-'LAMP5.MGE'

tab$broad[tab$cluster %in% c(18)]<-'INH'
tab$anno[tab$cluster %in% c(18)]<-'GABA'
tab$type[tab$cluster %in% c(18)]<-'LAMP5.CGE'

tab$broad[tab$cluster %in% c(39)]<-'INH'
tab$anno[tab$cluster %in% c(39)]<-'GABA'
tab$type[tab$cluster %in% c(39)]<-'CXCL14'

tab$broad[tab$cluster %in% c(40)]<-'INH'
tab$anno[tab$cluster %in% c(40)]<-'GABA'
tab$type[tab$cluster %in% c(40)]<-'HTR3A'

tab$broad[tab$cluster %in% c(29)]<-'INH'
tab$anno[tab$cluster %in% c(29)]<-'GABA'
tab$type[tab$cluster %in% c(29)]<-'VIP'

######NNC
tab$broad[tab$cluster %in% c(1,9)]<-'Glia'
tab$anno[tab$cluster %in% c(1,9)]<-'Oligo'
tab$type[tab$cluster %in% c(1,9)]<-paste0('Oligo.',c(1:2))

tab$broad[tab$cluster %in% c(8,20,27)]<-'Glia'
tab$anno[tab$cluster %in% c(8,20,27)]<-'Astro'
tab$type[tab$cluster %in% c(8,20,27)]<-paste0('Astro.',c(1:3))

tab$broad[tab$cluster %in% c(24)]<-'Glia'
tab$anno[tab$cluster %in% c(24)]<-'OPC'
tab$type[tab$cluster %in% c(24)]<-'OPC'

tab$broad[tab$cluster %in% c(44)]<-'Glia'
tab$anno[tab$cluster %in% c(44)]<-'OPC'
tab$type[tab$cluster %in% c(44)]<-'COP'

tab$broad[tab$cluster %in% c(2,6)]<-'Immune'
tab$anno[tab$cluster %in% c(2,6)]<-'Immune'
tab$type[tab$cluster %in% c(2,6)]<-paste0('Micro.',c(1:2))

tab$broad[tab$cluster %in% c(42)]<-'Immune'
tab$anno[tab$cluster %in% c(42)]<-'Immune'
tab$type[tab$cluster %in% c(42)]<-'Macro/Tcell'

tab$broad[tab$cluster %in% c(48)]<-'CSF'
tab$anno[tab$cluster %in% c(48)]<-'Ependy'
tab$type[tab$cluster %in% c(48)]<-'Ependy'

tab$broad[tab$cluster %in% c(34,36,62)]<-'CSF'
tab$anno[tab$cluster %in% c(34,36,62)]<-'Choroid'
tab$type[tab$cluster %in% c(34,36,62)]<-paste0('CP.',c(1:3))

tab$broad[tab$cluster %in% c(19,47)]<-'Vascular'
tab$anno[tab$cluster %in% c(19,47)]<-'Vascular'
tab$type[tab$cluster %in% c(19,47)]<-paste0('Endo.',c(1:2))

tab$broad[tab$cluster %in% c(30)]<-'Vascular'
tab$anno[tab$cluster %in% c(30)]<-'Vascular'
tab$type[tab$cluster %in% c(30)]<-'PC/SMC'

tab$broad[tab$cluster %in% c(41)]<-'Vascular'
tab$anno[tab$cluster %in% c(41)]<-'Vascular'
tab$type[tab$cluster %in% c(41)]<-'VLMC'


sce$broad.type<-factor(tab$broad[match(sce$k_5_louvain,tab$cluster)],
                       levels=c("EXC", "INH", "Glia", "Immune",
                                "CSF", "Vascular"))
sce$cell.type<-factor(tab$anno[match(sce$k_5_louvain,tab$cluster)],
                       levels=c("GC", "CA2-4", "CA1", "ProS/Sub", "L2/3", "L5",
                                "L6/6b", "HATA/AHi", "Thal", "Cajal", "GABA",
                                "Oligo", "Astro", "OPC", "Immune", "Ependy",
                                "Choroid", "Vascular"))
sce$fine.type<-factor(tab$type[match(sce$k_5_louvain,tab$cluster)],
                      levels=c("GC.1", "GC.2", "GC.3", "GC.4", "GC.5", "MC",
                               "CA3.1", "CA3.2", "CA2", "CA1", "ProS", "Sub.1",
                               "Sub.2", "L2/3.1", "L2/3.2", "L2/3.3", "L2/3.4",
                               "L2/3.5", "L2/3.6", "L5.1", "L5.2", "L6.1", "L6.2",
                               "L6b", "HATA", "AHi.1", "AHi.2", "AHi.3", "AHi.4",
                                "Thal","Cajal", "PENK", "SST", "CORT", "PV.FS",
                               "CRABP1", "C1QL1", "LAMP5.MGE", "LAMP5.CGE",
                               "CXCL14", "HTR3A", "VIP", "Oligo.1", "Oligo.2",
                               "Astro.1", "Astro.2", "Astro.3", "OPC", "COP",
                               "Micro.1", "Micro.2", "Macro/Tcell", "Ependy",
                               "CP.1", "CP.2", "CP.3", "Endo.1", "Endo.2",
                               "PC/SMC", "VLMC"))

features<-c('PROX1','ABI3BP','TSPAN18','ST18','RGS14','FNDC1','CHRM5','GFRA1','COL12A1',
            'TLE4','THEMIS','RORB','SATB2','CUX2','ESR1','MOXD1','DCSTAMP','NPFFR2','PAPPA2',
            'SHOX2','NDNF','PENK','SST','CORT','BTBD11','KCNAB3','C1QL1','CRABP1',
            'LAMP5','KIT','CXCL14','CNR1','VIP',
            'ETNPPL','MOBP','VCAN','GPR17','C3','SKAP1','CFAP73','PRLR',
            'VWF','ACTA2','CEMIP')

features<-c('PROX1','TSPAN18','FNDC1','CHRM5','GFRA1','COL12A1',
            'TLE4','THEMIS','RORB','SATB2','CUX2','ESR1','MOXD1','DCSTAMP','NPFFR2','PAPPA2',
            'SHOX2','NDNF','PENK','SST','CORT','BTBD11','KCNAB3','C1QL1','CRABP1',
            'LAMP5','KIT','CXCL14','CNR1','VIP',
            'ETNPPL','MOBP','PDGFRA','GPR17','C3','SKAP1','CFAP73','PRLR',
            'VWF','ACTA2','CEMIP')


plotGroupedHeatmap(sce,features=features,group='fine.type',cluster_cols=F,cluster_rows=F,color = brewer.pal(n=9,'Purples'),zlim=c(0,5))

levels<-c(1:16,22,18:21,17,23,25,24,27,26,28:52)
sce$cell.type2<-factor(sce$cell.type,levels=levels(sce$cellType)[levels])
levels(sce$cell.type2)[44]<-'OPC'
heat<-sce[features,]


heat <- summarizeAssayByGroup(
    heat,
    ids=heat$fine.type,
    subset.row = NULL,
    subset.col = NULL,
    statistics = c("mean", "sum", "num.detected", "prop.detected", "median"),
    store.number = "ncells",
    threshold = 0,
    BPPARAM = SerialParam(),
    assay.type='logcounts'
)



pdf(file=here::here('plots','figures','figure_3','heatmap_marker genes_legend_presentation.pdf'),w=18,h=12)
pheatmap(as.data.frame(assay(heat,'mean')),cluster_cols=F,cluster_rows=F,scale='none',
         color=brewer.pal(9,'RdPu'),breaks=seq(0,5,by=5/9),
         annotation_col=annotation,
         #annotation_row=annotation3,
         annotation_colors = palettes,annotation_legend=F,legend=T)
dev.off()

pdf(file=here::here('plots','figures','figure_3','heatmap_markers_legend.pdf'),w=5,h=2.5)
pheatmap(logcounts(heat),cluster_cols=F,cluster_rows=F,scale='column',
         color=magma(25),breaks=seq(-4,4,by=0.33),
         annotation_col=as.data.frame(annotation),
         annotation_colors = palettes,annotation_legend=F)
dev.off()

features<-c('PROX1','ABI3BP','TSPAN18','ST18','RGS14','FNDC1','CHRM5','GFRA1','COL12A1',
            'TLE4','THEMIS','RORB','SATB2','CUX2','ESR1','MOXD1','DCSTAMP','NPFFR2','PAPPA2',
            'SHOX2','NDNF','PENK','SST','CORT','BTBD11','KCNAB3','C1QL1','CRABP1',
            'LAMP5','KIT','CXCL14','CNR1','VIP',
            'ETNPPL','MOBP','VCAN','GPR17','C3','SKAP1','CFAP73','PRLR',
            'VWF','ACTA2','CEMIP')

features<-c('PROX1','TSPAN18','FNDC1','CHRM5','GFRA1','COL12A1',
            'TLE4','THEMIS','RORB','SATB2','CUX2','ESR1','MOXD1','DCSTAMP','NPFFR2','PAPPA2',
            'SHOX2','NDNF','PENK','SST','CORT','BTBD11','KCNAB3','C1QL1','CRABP1',
            'LAMP5','KIT','CXCL14','CNR1','VIP',
            'ETNPPL','MOBP','PDGFRA','GPR17','C3','SKAP1','CFAP73','PRLR',
            'VWF','ACTA2','CEMIP')


plotGroupedHeatmap(sce,features=features,group='fine.type',cluster_cols=F,cluster_rows=F,color = brewer.pal(n=9,'Purples'),zlim=c(0,5))

levels<-c(1:16,22,18:21,17,23,25,24,27,26,28:52)
sce$cell.type2<-factor(sce$cell.type,levels=levels(sce$cellType)[levels])
levels(sce$cell.type2)[44]<-'OPC'
heat<-sce[features,]


heat <- summarizeAssayByGroup(
    heat,
    ids=heat$fine.type,
    subset.row = NULL,
    subset.col = NULL,
    statistics = c("mean", "sum", "num.detected", "prop.detected", "median"),
    store.number = "ncells",
    threshold = 0,
    BPPARAM = SerialParam(),
    assay.type='logcounts'
)



pdf(file=here::here('plots','figures','figure_3','heatmap_marker genes_legend_presentation.pdf'),w=18,h=12)
pheatmap(as.data.frame(assay(heat,'mean')),cluster_cols=F,cluster_rows=F,scale='none',
         color=brewer.pal(9,'RdPu'),breaks=seq(0,5,by=5/9),
         annotation_col=annotation,
         #annotation_row=annotation3,
         annotation_colors = palettes,annotation_legend=F,legend=T)
dev.off()

pdf(file=here::here('plots','figures','figure_3','heatmap_markers_legend.pdf'),w=5,h=2.5)
pheatmap(logcounts(heat),cluster_cols=F,cluster_rows=F,scale='column',
         color=magma(25),breaks=seq(-4,4,by=0.33),
         annotation_col=as.data.frame(annotation),
         annotation_colors = palettes,annotation_legend=F)
dev.off()

# Obtain the column metadata from SpatialExperiment
df <- as.data.frame(colData(sce))

# Obtain the column metadata from SpatialExperiment
df <- as.data.frame(colData(sce))

# Create a data frame with counts of each cluster for each sample (brnum)
count_df <- df %>%
    count(cell.type, brnum) %>%
    group_by(brnum) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
# Set the factor levels of sample_id based on the order in neuron_order
count_df$brnum <- factor(count_df$brnum, levels = neuron_order)

palette<-c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
           "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
           "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
           "#C7C7C7", "#BCBD22", "#DBDB8D")


neuron_order<-c("Br2720", "Br6522", "Br8492", "Br8325", "Br3942",
                "Br8667", "Br6471", "Br6432", "Br6423", "Br2743")


# Aggregate data by broad
agg_df <- count_df %>%
    group_by(cell.type) %>%
    summarise(n = sum(n), proportion = sum(n) / sum(count_df$n))

# Generate plot
pdf(file=here::here('plots','figures','figure_3','barplot_proportion_total_celltype.pdf'),
    h=2.5,w=4)
ggplot(agg_df, aes(x = "Total", y = proportion, fill = cell.type)) +
    geom_bar(stat = 'identity') +
    labs(x = "", y = 'Overall Proportion', fill = 'broad') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(colour='black')) +
    scale_fill_manual(values = palettes[[1]])
dev.off()

pdf(file=here::here('plots','figures','figure_3','barplot_proportion_sort.pdf'),
    w=2.5,h=4)
ggplot(count_df, aes(x = sort, y = proportion, fill = cell.type)) +
    geom_bar(stat = 'identity', position = 'fill') +
    labs(x = 'sort', y = 'proportion', fill = 'cell type') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(colour='black'),
          legend.position='none') +
    scale_fill_manual(values=palettes[[1]])
dev.off()


pdf(file=here::here('plots','figures','figure_3','barplot_proportion_celltype.pdf'),
    w=2.5,h=4)
ggplot(count_df, aes(x = brnum, y = proportion, fill = cell.type)) +
    geom_bar(stat = 'identity', position = 'fill') +
    labs(x = 'donor', y = 'proportion', fill = 'cell type') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(colour='black'),
          legend.position='none') +
    scale_fill_manual(values=palettes[[1]])
dev.off()

# Create a data frame with counts of each cluster for each sample (sort)
broad_df <- df %>%
    count(broad.type, sort, brnum) %>%
    group_by(sort) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()


pdf(file=here::here('plots','figures','figure_3','barplot_proportion_donor_broad.pdf'),
    w=2.5,h=4)
ggplot(broad_df, aes(x = brnum, y = proportion, fill = broad.type)) +
    geom_bar(stat = 'identity', position = 'fill') +
    labs(x = 'donor', y = 'proportion', fill = 'broad type') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(colour='black'),
          legend.position='none') +
    scale_fill_manual(values=palettes[[2]])
dev.off()

pdf(file=here::here('plots','figures','figure_3','barplot_proportion_sort_broad.pdf'),
    w=2.5,h=4)
ggplot(broad_df, aes(x = sort, y = proportion, fill = broad.type)) +
    geom_bar(stat = 'identity', position = 'fill') +
    labs(x = 'sort', y = 'proportion', fill = 'broad type') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(colour='black'),
          legend.position='none') +
    scale_fill_manual(values=palettes[[2]])
dev.off()


combined_df$facet <- ifelse(combined_df$sort == "Overall", "Overall", "Individual")

# Create the barplot using ggplot2
pdf(file=here::here('plots','figures','figure_3','barplot_proportion_broad.pdf'),
    h=3.5,w=2)
ggplot(broad_df, aes(x = sort, y = proportion, fill = broad)) +
    geom_bar(stat = 'identity', position = 'fill') +
    labs(x = 'Sort', y = 'Proportion', fill = 'broad') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(colour='black'),
          legend.position='none') +
    geom_text(data = unique(select(broad_df, sort, broad)),
              aes(x = sort, y = 0, label = sort),
              angle = 0, hjust = 0.5, vjust = 1, check_overlap = TRUE)+
    scale_fill_manual(values=palette2)


dev.off()




# pdf(file=here::here('plots','figures','figure_3','cluster_UMAP.pdf'),h=3.5,w=3.75)
# plot<-plotUMAP(sce,text_by='cell.type',colour_by='cell.type',point_size=0.001,add_legend=F)+
#     theme(axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank(),
#           text = element_text(colour='black',size=11))
# rasterize(plot,dpi=1000)
# dev.off()



#
# # Load necessary packages
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(Polychrome)
#
# # Obtain the column metadata from SpatialExperiment
# df <- as.data.frame(colData(sce))
#
# # Create a data frame with counts of each cluster for each sample (brnum)
# count_df <- df %>%
#     count(annotation, brnum) %>%
#     group_by(brnum) %>%
#     mutate(proportion = n / sum(n)) %>%
#     ungroup()
#
# neuron_order<-c("Br2720", "Br6522", "Br8492", "Br8325", "Br3942",
#                 "Br8667", "Br6471", "Br6432", "Br6423", "Br2743")
#
# # Set the factor levels of sample_id based on the order in neuron_order
# count_df$brnum <- factor(count_df$brnum, levels = neuron_order)
#
# palette<-c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
#                        "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
#                        "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
#                        "#C7C7C7", "#BCBD22", "#DBDB8D")
#
# combined_df$facet <- ifelse(combined_df$brnum == "Overall", "Overall", "Individual")
#
# # Create the barplot using ggplot2
# pdf(file=here::here('plots','figures','figure_3','barplot_proportion.pdf'),
#     h=3.5,w=3.5)
# ggplot(count_df, aes(x = brnum, y = proportion, fill = annotation)) +
#     geom_bar(stat = 'identity', position = 'fill') +
#     labs(x = 'Donor', y = 'Proportion', fill = 'annotation') +
#     theme_minimal() +
#     theme(axis.text.x = element_blank(),
#           text = element_text(colour='black'),
#           legend.position='none') +
#     geom_text(data = unique(select(count_df, brnum, annotation)),
#               aes(x = brnum, y = 0, label = brnum),
#               angle = 0, hjust = 0.5, vjust = 1, check_overlap = TRUE)+
#     scale_fill_manual(values=palette)
# dev.off()
#
# sce$cellType_collapsed<-tab$annoType[match(sce$k_5_louvain_1.2,tab$cluster)]
# features=c('C3','CSF1R')
# plotDots(sce,group='k_5_louvain_initial',features=rownames(marks[[45]][1:25,]),color=c('white','black')) #+
#     #scale_y_discrete(limits=rev(features)) +
#     scale_x_discrete(limits=rev(levels(spe$cluster)))+
#     theme(axis.text.x = element_text(angle = 90,vjust=0.75),
#           text = element_text(size = 14))+ coord_flip()+
#     labs(x= "Cell type", y = "Gene")
#
# levels(sce$cell.type)
#
#        "CP1"         "CP2"         "Endo"
#      "PC"                "SMC"                "VLMC" "VIP"
#
#      # [1] "GC.1"        "GC.2"        "GC.3"        "GC.4"        "GC.5"        "MC"
#      # [7] "CA3.1"       "CA3.2"       "CA2"         "CA1"         "ProS"        "Sub.1"
#      # [13] "Sub.2"       "L6.2"        "L6.1"        "L6b"         "L5.2"        "L5.1"
#      # [19] "L2/3.1"      "L2/3.2"      "L2/3.3"      "L2/3.4"      "L2/3.5"      "L2/3.6"
#      # [25] "HATA"        "AHi.1"       "AHi.2"       "AHi.3"       "AHi.4"       "Thal"
#      # [31] "Cajal"       "PENK"        "SST"         "CORT"        "PV.FS"       "C1QL1"
#      # [37] "CRABP1"      "LAMP5.MGE"   "LAMP5.CGE"   "CXCL14"      "HTR3A"       "VIP"
#      # [43] "Astro.1"     "Astro.2"     "Astro.3"     "Oligo.1"     "Oligo.2"     "OPC"
#      # [49] "COP"         "Micro.1"     "Micro.2"     "Macro/Tcell" "Ependy"      "CP.1"
#      # [55] "CP.2"        "CP.3"        "Endo.2"      "Endo.1"      "PC/SMC"      "VLMC"
#
#  sce$cell.type<-factor(sce$cell.type,levels=levels)
#
#  features<-c('SYT1','SLC17A7','PROX1','TSPAN18','FNDC1','ST8SIA2','RPRM','SATB2'
#              'TLE4','RORB','CUX2','ESR1','NPFFR2','MOXD1','PAPPA2','DCSTAMP',
#              'SHOX2','NDNF','GAD2','ETNPPL','MOBP','VCAN','C3','CFAP73','PRLR','EBF1')
#  #features<-rownames(marks.anno[[18]])[1:40]
#  pdf(file=here::here('plots','figures','figure_3','dotplot.pdf'),h=7,w=16)
#
#  features<-c('SERPINA5','RYBP','SLC38A2','FEM1B','PYDC1','APP','PSEN1','PSEN2','APOE',
#              'ADAM10','ANKMY2')#,'ATP5F1')
#
#  ADAM10, ANKMY2, ATP5F1
#  plotDots(sce,group='cell.type',features=features,color=c('white','black')) +
#      scale_y_discrete(limits=rev(features)) +
#      #scale_x_discrete(limits=rev(levels(sce$cell.type)))+
#      theme(axis.text.x = element_text(angle = 90,vjust=0.75),
#            text = element_text(size = 14))+
#      labs(x= "cell type", y = "gene")
#  dev.off()
#
#  50,44,45,41,49
