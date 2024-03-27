setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("scuttle"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("scry"))
suppressPackageStartupMessages(library("jaffelab"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("BiocSingular"))

###load data
load(here::here("snRNAseq_hpc",'processed-data','sce','sce_clustered_round2.rda'), verbose = TRUE)
sce$cluster<-sce$k_5_louvain_initial


sce$neuron<-ifelse(sce$k_5_louvain_initial %in% c(1,2,9,10,20,21,26,28,
                                                  31,34,37,38,40,42,43,
                                                  59),T,F)
sce$discard_round_3<-ifelse(sce$k_5_louvain_initial %in% c(6,8,46,37),T,F)

table(sce$discard)
# FALSE  TRUE
# 3926 76668

###first a UMAP by cluster
pdf(file=here::here('plots','figures','supp_figures',
                    'clusteringRound2_UMAP_byCluster.pdf'),h=3,w=3)
p<-plotUMAP(sce,text_by='cluster',colour_by='cluster',point_size=0.01,add_legend=F,text_size=3.75)
rasterize(p,dpi=500,layers='Point')
dev.off()

###now boxplot of SYT1 expression by cluster2
pdf(file=here::here('plots','figures','supp_figures',
                    'clusteringRound2_boxplot_SYT1_byCluster.pdf'),h=3,w=5)
ggcells(sce, mapping=aes(x=cluster, y=SYT1,fill=neuron)) +
    geom_boxplot(outlier.size = 0.5)+theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 11,colour='black'))+
    theme(legend.position='bottom',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_fill_manual(values=c('green','blue'))
dev.off()

###now a UMAP by cluster coloured by neurons
pdf(file=here::here('plots','figures','supp_figures',
                    'clusteringRound2_UMAP_byCluster_neurons.pdf'),h=3,w=4)
p<-plotUMAP(sce,text_by='cluster',colour_by='neuron',
            point_size=0.01,add_legend=T,text_size=3.75)+
                scale_color_manual(values=c('green','blue'))+labs(color='neuron')
rasterize(p,dpi=500,layers='Point')
dev.off()



###now boxplot of detected genes by cluster
pdf(file=here::here('plots','figures','supp_figures',
                    'clusteringRound2_boxplot_detectedGenes_byCluster.pdf'),h=3,w=5)
ggcells(sce, mapping=aes(x=cluster, y=detected,fill=neuron)) +
    geom_boxplot(outlier.size = 0.5)+theme(axis.text.x = element_text(angle = 90),
                                           text=element_text(size = 11,colour='black'))+
    theme(legend.position='bottom',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    scale_fill_manual(values=c('green','blue'))+
    ylab('detected genes')+scale_y_log10()

dev.off()

##roundabout but easy way of getting ordered clusters for doublet plot
p1<-plotGroupedHeatmap(sce,group = 'cluster',
                       features=c('C3','CSF1R',

                           'AQP4','GFAP','ETNPPL','MOBP',
                                  'MBP','MOG','PDGFRA','VCAN'),
                       center=T,scale=T)
lab_order<-as.numeric(p1$tree_col$order)
sce$cluster_ordered<-factor(sce$cluster,levels=lab_order)

##doublet dotplot
pdf(file=here::here('plots','figures','supp_figures',
                    'clusteringRound2_dotPlot_markerGenes_doublets.pdf'),h=6,w=17)
plotDots(sce,group='cluster_ordered',
         features=c('PDGFRA',
                    'AQP4','ETNPPL',
                    'MOBP','MOG',
                    'C3','CSF1R'),
         color=c('white','black')) +
    #scale_y_discrete(limits=rev(features)) +
    #scale_x_discrete(limits=rev(levels(spe$cluster)))+
    theme(axis.text.x = element_text(angle = 90,vjust=0.75),
          text = element_text(size = 20))+
    labs(x= "Cell type", y = "Gene",size='proportion',color='mean')
dev.off()

##discard UMAP
pdf(file=here::here('plots','figures','supp_figures',
                    'clusteringRound2_UMAP_byCluster_discardRound3.pdf'),h=3,w=4.5)
p<-plotUMAP(sce,text_by='cluster',colour_by='discard_round_3',
            point_size=0.01,add_legend=T,text_size=3.75)+
    labs(color='discard_round_3')+scale_color_manual(values=c('#859ECA','#FF9E4A'))
rasterize(p,dpi=500,layers='Point')
dev.off()

##save data 
save(sce,file=here::here("snRNAseq_hpc",'processed-data','sce','sce_clustered_round2.rda'))

##remove sce
rm(sce)

##load next round sce
load(here::here("snRNAseq_hpc",'processed-data','sce','sce_clustered_round3.rda'), verbose = TRUE)


# ##discard clusters, make new embedding and UMAP colored by discard round 3, but numbered by new cluster round
#pdf(file=here::here('plots','figures','supp_figures',
#                     'clusteringRound3_UMAP_byCluster_discardRound3.pdf'),h=3,w=4.5)
# p<-plotUMAP(sce,colour_by='discard_round_3',text_by='k_5_louvain',
#             point_size=0.01,add_legend=T,text_size=3.75)+
#     labs(color='discard_round_3')+scale_color_manual(values='#859ECA')
# rasterize(p,dpi=500,layers='Point')
# dev.off()


sce$neuron<-ifelse(sce$k_5_louvain %in% c(1,2,6,8,9,19,20,23,24,27,30,34,36,41,42,44,47,48,62),F,T)
sce$cluster<-sce$k_5_louvain

###now boxplot of SYT1 expression by cluster3
pdf(file=here::here('plots','figures','supp_figures',
                    'clusteringRound3_boxplot_SYT1_byCluster.pdf'),h=3,w=5)
ggcells(sce, mapping=aes(x=k_5_louvain, y=SYT1,fill=neuron)) +
    geom_boxplot(outlier.size = 0.5)+theme(axis.text.x = element_text(angle = 90),
                                           text=element_text(size = 11,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+xlab('cluster3')+
    scale_fill_manual(values=c('blue','green'))
dev.off()

###now a UMAP by cluster coloured by neurons
pdf(file=here::here('plots','figures','supp_figures',
                    'clusteringRound3_UMAP_byCluster_neurons.pdf'),h=3,w=4)
p<-plotUMAP(sce,text_by='cluster',colour_by='neuron',
            point_size=0.01,add_legend=T,text_size=3.75)+
    scale_color_manual(values=c('blue','green'))+labs(color='neuron')
rasterize(p,dpi=500,layers='Point')
dev.off()



###now boxplot of detected genes by cluster
pdf(file=here::here('plots','figures','supp_figures',
                    'clusteringRound3_boxplot_detectedGenes_byCluster.pdf'),h=3,w=6)
ggcells(sce, mapping=aes(x=cluster, y=detected,fill=neuron)) +
    geom_boxplot(outlier.size = 0.5)+theme(axis.text.x = element_text(angle = 90),
                                           text=element_text(size = 11,colour='black'))+
    theme(legend.position='bottom',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    scale_fill_manual(values=c('blue','green'))+
    ylab('detected genes')+scale_y_log10()

dev.off()


##roundabout but easy way of getting ordered clusters for doublet plot
p1<-plotGroupedHeatmap(sce,group = 'cluster',
                       features=c('C3',
                                  'AQP4','GFAP','ETNPPL','MOBP',
                                  'MBP','MOG','PDGFRA','VCAN','GPR17'),
                       center=T,scale=T)
lab_order<-as.numeric(p1$tree_col$order)
sce$cluster_ordered<-factor(sce$cluster,levels=lab_order)

##doublet dotplot
features=c('MOBP','MOG',
           'AQP4','ETNPPL',
           'PDGFRA','VCAN','GPR17')
pdf(file=here::here('plots','figures','supp_figures',
                    'clusteringRound3_dotPlot_markerGenes_doublets.pdf'),h=6,w=17)
plotDots(sce,group='cluster_ordered',
         features=features,
         color=c('white','black')) +
    scale_y_discrete(limits=rev(features)) +
    #scale_x_discrete(limits=rev(levels(spe$cluster)))+
    theme(axis.text.x = element_text(angle = 90,vjust=0.75),
          text = element_text(size = 20))+
    labs(x= "Cell type", y = "Gene",size='proportion',color='mean')
dev.off()

sce$discard_round_4<-ifelse(sce$k_5_louvain %in% c(17,23),T,F)
pdf(file=here::here('plots','figures','supp_figures',
                    'clusteringRound2_UMAP_byCluster_discardRound4.pdf'),h=3,w=4.5)
p<-plotUMAP(sce,text_by='k_5_louvain',colour_by='discard_round_4',
            point_size=0.01,add_legend=T,text_size=3.75)+
    labs(color='discard_round_4')+scale_color_manual(values=c('#859ECA','#FF9E4A'))
rasterize(p,dpi=500,layers='Point')
dev.off()

##save data
save(sce,file=here::here("snRNAseq_hpc",'processed-data','sce','sce_clustered_round3.rda'))

##remove sce
rm(sce)

##load next round data
load(file=here::here("snRNAseq_hpc",'processed-data','sce','sce_final.rda'))


##discard clusters, make new embedding and UMAP colored by discard round 3, but numbered by new cluster round
pdf(file=here::here('plots','figures','supp_figures',
                    'clusteringRound4_UMAP_byCluster_discardRound4.pdf'),h=3,w=5.5)
p<-plotUMAP(sce,colour_by='cell.class',text_by='cell.class',
            point_size=0.01,add_legend=T,text_size=3.75)+
    labs(color='fine.cell.class')+scale_color_manual(values=sn.fine.palette)
rasterize(p,dpi=500,layers='Point')
dev.off()
