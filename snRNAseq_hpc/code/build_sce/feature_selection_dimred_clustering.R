library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("here")
library("sessioninfo")
library("harmony")
library("BiocSingular")


## Load sce
load(here("snRNAseq_hpc","processed-data", "sce", "sce_post_qc.rda"))
dim(sce)

#Remove nuclei from problematic sample
sce$discard<-ifelse(sce$Sample %in% "17c-scp" & sce$detected < 5000,T,F)
sce<-sce[,sce$discard==F]

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

abline(v = 3000, lty = 2, col = "blue")
dev.off()

hdg<-rownames(counts(sce))[1:3000]
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
res <- scater::runPCA(res,
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
mnn<-batchelor::reducedMNN(reducedDim(sce,'PCA'),batch=sce$round,k=50)

reducedDim(sce,'MNN')<-mnn$corrected
set.seed(2351)
message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "MNN")

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
clust_5 <- igraph::cluster_louvain(snn.gr,resolution=1.2)$membership
table(clust_5)
sce$k_5_louvain_initial<-factor(clust_5)

#pdf(here("snRNAseq_hpc","plots","UMAPs_corrected.pdf"))
plotUMAP(sce,text_by='k_5_louvain_initial',
         colour_by='k_5_louvain_initial',point_size=0.25)


##flag some low-quality neuron clusters
sce$discard<-ifelse(sce$k_5_louvain_initial %in% c(6,8,18,25,35,39),T,F)
sce<-sce[,sce$discard==F]

##make neuron vs non-neuron comparisons
index<-c(1,27,21,9,2,41,33,49,38,57,29,20)
sce$neuron<-ifelse(sce$k_5_louvain_initial %in% index,F,T)


##Rerun all steps
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

abline(v = 3000, lty = 2, col = "blue")
dev.off()

hdg<-rownames(counts(sce))[1:3000]
res<-sce[rownames(counts(sce)) %in% hdg,]
set.seed(913)
res <- nullResiduals(res,
                     fam = "poisson",
                     type = "pearson")

set.seed(915)
message("running PCA - ", Sys.time())
res <- scater::runPCA(res,
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
mnn<-batchelor::reducedMNN(reducedDim(sce,'PCA'),batch=sce$round,k=50)

reducedDim(sce,'MNN')<-mnn$corrected
set.seed(1122)
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
sce$k_5_louvain<-factor(clust_5)


message("running buildSNNGraph - ", Sys.time())
snn.gr <- buildSNNGraph(sce, k = 10, use.dimred = "MNN",type='jaccard')

set.seed(100)
clust_10 <- igraph::cluster_louvain(snn.gr)$membership
table(clust_10)
sce$k_10_louvain<-factor(clust_10)

set.seed(100)
message("running buildSNNGraph - ", Sys.time())
snn.gr <- buildSNNGraph(sce, k = 5, use.dimred = "MNN",type='jaccard')
clust_5 <- igraph::cluster_louvain(snn.gr,resolution=1.75)$membership
table(clust_5)
sce$k_5_louvain<-factor(clust_5)

marks<-findMarkers(sce,groups=sce$k_5_louvain,pval.type='all',direction='up',test.type='wilcox')





tab<-data.frame('cluster'=1:52,'broad'=rep(NA,52),
                'anno'=rep(NA,52),
                'type'=rep(NA,52))
########Neurons
##Granule cells
tab$broad[tab$cluster %in% c(10,11,29,37,42)]<-'EXC'
tab$anno[tab$cluster %in% c(10,11,29,37,42)]<-'GC'
tab$type[tab$cluster %in% c(10,11,29,37,42)]<-paste0('GC.',c(1:5))

##Mossy cells
tab$broad[tab$cluster %in% c(18)]<-'EXC'
tab$anno[tab$cluster %in% c(18)]<-'CA2-4'
tab$type[tab$cluster %in% c(18)]<-'MC'
##CA3
tab$broad[tab$cluster %in% c(9,31)]<-'EXC'
tab$anno[tab$cluster %in% c(9,31)]<-'CA2-4'
tab$type[tab$cluster %in% c(9,31)]<-paste0('CA3.',c(1:2))
##CA2
tab$broad[tab$cluster %in% c(34)]<-'EXC'
tab$anno[tab$cluster %in% c(34)]<-'CA2-4'
tab$type[tab$cluster %in% c(34)]<-'CA2'

##CA1
tab$broad[tab$cluster %in% c(13)]<-'EXC'
tab$anno[tab$cluster %in% c(13)]<-'CA1'
tab$type[tab$cluster %in% c(13)]<-'CA1'

##Sub_sup
tab$broad[tab$cluster %in% c(23)]<-'EXC'
tab$anno[tab$cluster %in% c(23)]<-'Sub'
tab$type[tab$cluster %in% c(23)]<-'Sub.1'
##Sub_deep
tab$broad[tab$cluster %in% c(19)]<-'EXC'
tab$anno[tab$cluster %in% c(19)]<-'Sub'
tab$type[tab$cluster %in% c(19)]<-'Sub.2'

##RHP
##L2/3
tab$broad[tab$cluster %in% c(5,25,40,36,38,12)]<-'EXC'
tab$anno[tab$cluster %in% c(5,25,40,36,38,12)]<-'L2/3'
tab$type[tab$cluster %in% c(5,25,40,36,38,12)]<-paste0('L2/3.',c(1:6))

##L5
tab$broad[tab$cluster %in% c(21)]<-'EXC'
tab$anno[tab$cluster %in% c(21)]<-'L5'
tab$type[tab$cluster %in% c(21)]<-'L5'

##L6
tab$broad[tab$cluster %in% c(14,28)]<-'EXC'
tab$anno[tab$cluster %in% c(14,28)]<-'L6/6b'
tab$type[tab$cluster %in% c(14,28)]<-paste0('L6.',c(1,2))

##L6b
tab$broad[tab$cluster %in% c(20)]<-'EXC'
tab$anno[tab$cluster %in% c(20)]<-'L6/6b'
tab$type[tab$cluster %in% c(20)]<-'L6b'

##Thalamus
tab$broad[tab$cluster %in% c(4)]<-'EXC'
tab$anno[tab$cluster %in% c(4)]<-'Thal'
tab$type[tab$cluster %in% c(4)]<-'Thal'



##AHi
tab$broad[tab$cluster %in% c(50,44,45,41,49)]<-'EXC'
tab$anno[tab$cluster %in% c(50,44,45,41,49)]<-'HATA/AHi'
tab$type[tab$cluster %in% c(50,44,45,41,49)]<-paste0('AHi.',c(1:5))
##HATA
tab$broad[tab$cluster %in% c(41)]<-'EXC'
tab$anno[tab$cluster %in% c(41)]<-'HATA/AHi'
tab$type[tab$cluster %in% c(41)]<-'HATA'

##Cajal Retzius
tab$broad[tab$cluster %in% c(47)]<-'EXC'
tab$anno[tab$cluster %in% c(47)]<-'Cajal'
tab$type[tab$cluster %in% c(47)]<-'Cajal'

##GABAergics
tab$broad[tab$cluster %in% c(48)]<-'INH'
tab$anno[tab$cluster %in% c(48)]<-'GABA'
tab$type[tab$cluster %in% c(48)]<-'MEIS2'

tab$broad[tab$cluster %in% c(15)]<-'INH'
tab$anno[tab$cluster %in% c(15)]<-'GABA'
tab$type[tab$cluster %in% c(15)]<-'PV'

tab$broad[tab$cluster %in% c(46)]<-'INH'
tab$anno[tab$cluster %in% c(46)]<-'GABA'
tab$type[tab$cluster %in% c(46)]<-'CRABP1'

tab$broad[tab$cluster %in% c(8)]<-'INH'
tab$anno[tab$cluster %in% c(8)]<-'GABA'
tab$type[tab$cluster %in% c(8)]<-'SST'

tab$broad[tab$cluster %in% c(39)]<-'INH'
tab$anno[tab$cluster %in% c(39)]<-'GABA'
tab$type[tab$cluster %in% c(39)]<-'CORT'

tab$broad[tab$cluster %in% c(3)]<-'INH'
tab$anno[tab$cluster %in% c(3)]<-'GABA'
tab$type[tab$cluster %in% c(3)]<-'LAMP5.CHST9'

tab$broad[tab$cluster %in% c(16)]<-'INH'
tab$anno[tab$cluster %in% c(16)]<-'GABA'
tab$type[tab$cluster %in% c(16)]<-'LAMP5.KIT'

tab$broad[tab$cluster %in% c(32)]<-'INH'
tab$anno[tab$cluster %in% c(32)]<-'GABA'
tab$type[tab$cluster %in% c(32)]<-'HTR3A'

tab$broad[tab$cluster %in% c(26)]<-'INH'
tab$anno[tab$cluster %in% c(26)]<-'GABA'
tab$type[tab$cluster %in% c(26)]<-'SST'

######NNC
tab$broad[tab$cluster %in% c(1,7)]<-'NNC'
tab$anno[tab$cluster %in% c(1,7)]<-'Oligo'
tab$type[tab$cluster %in% c(1,7)]<-paste0('Oligo.',c(1:2))

tab$broad[tab$cluster %in% c(6,24)]<-'NNC'
tab$anno[tab$cluster %in% c(6,24)]<-'Astro'
tab$type[tab$cluster %in% c(6,24)]<-paste0('Astro',c(1:2))

tab$broad[tab$cluster %in% c(22)]<-'NNC'
tab$anno[tab$cluster %in% c(22)]<-'OPC'
tab$type[tab$cluster %in% c(22)]<-'OPC'

tab$broad[tab$cluster %in% c(52)]<-'NNC'
tab$anno[tab$cluster %in% c(52)]<-'OPC'
tab$type[tab$cluster %in% c(52)]<-'COP'

tab$broad[tab$cluster %in% c(35)]<-'NNC'
tab$anno[tab$cluster %in% c(35)]<-'Ependy'
tab$type[tab$cluster %in% c(35)]<-'Ependy'

tab$broad[tab$cluster %in% c(2)]<-'NNC'
tab$anno[tab$cluster %in% c(2)]<-'Micro'
tab$type[tab$cluster %in% c(2)]<-'Micro'

tab$broad[tab$cluster %in% c(30,43)]<-'NNC'
tab$anno[tab$cluster %in% c(30,43)]<-'Choroid'
tab$type[tab$cluster %in% c(30,43)]<-paste0('CP',c(1:2))

tab$broad[tab$cluster %in% c(17)]<-'NNC'
tab$anno[tab$cluster %in% c(17)]<-'Vascular'
tab$type[tab$cluster %in% c(17)]<-'Endo'

tab$broad[tab$cluster %in% c(27)]<-'NNC'
tab$anno[tab$cluster %in% c(27)]<-'Vascular'
tab$type[tab$cluster %in% c(27)]<-'PC'

tab$broad[tab$cluster %in% c(33)]<-'NNC'
tab$anno[tab$cluster %in% c(33)]<-'Vascular'
tab$type[tab$cluster %in% c(33)]<-'VLMC'

tab$broad[tab$cluster %in% c(51)]<-'NNC'
tab$anno[tab$cluster %in% c(51)]<-'Vascular'
tab$type[tab$cluster %in% c(51)]<-'SMC'

sce$broad.type<-factor(tab$broad[match(sce$k_5_louvain_initial,tab$cluster)])
sce$annotation<-factor(tab$anno[match(sce$k_5_louvain_initial,tab$cluster)],
                       levels=c('GC','CA2-4','CA1','Sub',
                                'L2/3','L5','L6/6b','HATA/AHi',
                                'Thal','Cajal','GABA',
                                'Astro','Oligo','OPC',
                                'Micro','Ependy','Choroid','Vascular'))
sce$cell.type<-factor(tab$type[match(sce$k_5_louvain_initial,tab$cluster)])

pdf(file=here::here('plots','figures','figure_3','cluster_UMAP.pdf'),h=5,w=6)
plot<-plotUMAP(sce,text_by='cell.type',colour_by='COL5A2',point_size=0.1)+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          text = element_text(colour='black'))
rasterize(plot,layers = 'Point',dpi=1000)
dev.off()


# Load necessary packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(Polychrome)

# Obtain the column metadata from SpatialExperiment
df <- as.data.frame(colData(sce))

# Create a data frame with counts of each cluster for each sample (brnum)
count_df <- df %>%
    count(annotation, brnum) %>%
    group_by(brnum) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()

neuron_order<-c("Br2720", "Br6522", "Br8492", "Br8325", "Br3942",
                "Br8667", "Br6471", "Br6432", "Br6423", "Br2743")

# Set the factor levels of sample_id based on the order in neuron_order
count_df$brnum <- factor(count_df$brnum, levels = neuron_order)

palette<-c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
                       "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                       "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
                       "#C7C7C7", "#BCBD22", "#DBDB8D")

combined_df$facet <- ifelse(combined_df$brnum == "Overall", "Overall", "Individual")

# Create the barplot using ggplot2
pdf(file=here::here('plots','figures','figure_3','barplot_proportion.pdf'),
    h=3.5,w=3.5)
ggplot(count_df, aes(x = brnum, y = proportion, fill = annotation)) +
    geom_bar(stat = 'identity', position = 'fill') +
    labs(x = 'Donor', y = 'Proportion', fill = 'annotation') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(colour='black'),
          legend.position='none') +
    geom_text(data = unique(select(count_df, brnum, annotation)),
              aes(x = brnum, y = 0, label = brnum),
              angle = 0, hjust = 0.5, vjust = 1, check_overlap = TRUE)+
    scale_fill_manual(values=palette)
dev.off()

sce$cellType_collapsed<-tab$annoType[match(sce$k_5_louvain_1.2,tab$cluster)]

plotDots(sce,group='cell.type',features=features,color=c('white','black')) +
    #scale_y_discrete(limits=rev(features)) +
    scale_x_discrete(limits=rev(levels(spe$cluster)))+
    theme(axis.text.x = element_text(angle = 90,vjust=0.75),
          text = element_text(size = 14))+ coord_flip()+
    labs(x= "Cell type", y = "Gene")

levels(sce$cell.type)

       "CP1"         "CP2"         "Endo"
     "PC"                "SMC"                "VLMC" "VIP"

 levels<-c("GC.1","GC.2","GC.3","GC.4","GC.5","MC","CA3.1","CA3.2","CA2","CA1",
           "Sub.1","Sub.2","L6b","L6.1","L6.2","L5","L2/3.2","L2/3.1", "L2/3.3",
           "L2/3.4","L2/3.5","L2/3.6","HATA","AHi.2","AHi.3","AHi.4","AHi.5","Thal",
           "Cajal","PV","CRABP1","SST","CORT","LAMP5.CHST9","LAMP5.KIT","HTR3A","VIP",
           "MEIS2","Astro1","Astro2","Oligo.1","Oligo.2","OPC","COP","Micro","Ependy",
           "CP1","CP2","Endo","PC","SMC","VLMC")

 sce$cell.type<-factor(sce$cell.type,levels=levels)

 features<-c('SYT1','SLC17A7','PROX1','TSPAN18','FNDC1','ST8SIA2','RPRM','SATB2'
             'TLE4','RORB','CUX2','ESR1','NPFFR2','MOXD1','PAPPA2','DCSTAMP',
             'SHOX2','NDNF','GAD2','ETNPPL','MOBP','VCAN','C3','CFAP73','PRLR','EBF1')
 #features<-rownames(marks.anno[[18]])[1:40]
 pdf(file=here::here('plots','figures','figure_3','dotplot.pdf'),h=7,w=16)

 features<-c('SERPINA5','RYBP','SLC38A2','FEM1B','PYDC1','APP','PSEN1','PSEN2','APOE',
             'ADAM10','ANKMY2')#,'ATP5F1')

 ADAM10, ANKMY2, ATP5F1
 plotDots(sce,group='cell.type',features=features,color=c('white','black')) +
     scale_y_discrete(limits=rev(features)) +
     #scale_x_discrete(limits=rev(levels(sce$cell.type)))+
     theme(axis.text.x = element_text(angle = 90,vjust=0.75),
           text = element_text(size = 14))+
     labs(x= "cell type", y = "gene")
 dev.off()

 50,44,45,41,49
