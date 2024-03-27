#cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("scuttle"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("scry"))
suppressPackageStartupMessages(library("jaffelab"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("BiocSingular"))

##load data
load(here("snRNAseq_hpc",'processed-data','sce','sce_post_qc.rda'), verbose = TRUE)

##make sure feature names unique
rownames(sce)<-uniquifyFeatureNames(rowData(sce)$gene_id,rowData(sce)$gene_name)

##barplot showing  nuclei by sample, colored by discard status
# Convert SingleCellExperiment to dataframe
df <- as.data.frame(colData(sce))

# Calculate proportions
df <- df %>%
    group_by(sample_ID, sort) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup() %>%
    group_by(sample_ID) %>%
    mutate(total = sum(proportion)) %>%
    mutate(proportion = proportion / total)

# Create barplot
pdf(file=here::here('plots','figures','supp_figures','barplot_postDiscard.pdf'),h=4,w=8)
ggplot(df, aes(x = sample_ID, y = count)) +
    geom_bar(stat = "identity", position = "stack",fill='cornflowerblue') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text=element_text(colour='black',size=16),
          axis.text.x = element_text(angle = 90),
          legend.position='bottom') +
    ylab("nuclei") +
    xlab("sample") +
    labs(fill = "discard")+facet_wrap(~sort,scales='free_x')+
    geom_hline(yintercept=6000,linetype=2,color='magenta')+
    geom_hline(yintercept=9000,linetype=2,color='black')
dev.off()


##Now feature select/dimred/cluster
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

set.seed(125)
message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "MNN")


message("normalizing counts - ", Sys.time())
set.seed(100)
clust <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clust, min.mean=0.1)
sce <- logNormCounts(sce)

##make violin showing SYT1 vs #det genes
pdf(file=here::here('plots','figures','supp_figures','SYT1_vs_detected_violin_snRNAseq.pdf'),w=4,h=3)
p1<-plotColData(sce, x = "Sample", y = "detected",
                colour_by = "SYT1",point_alpha=1,point_size=0.001) +
    scale_y_log10() + ylab('detected genes')+xlab('sample')+
    facet_wrap(~sce$sort,scales='free_x')+
    ggtitle("Total detected genes")+
    theme(axis.text.x = element_text(size=8, angle = 90),
          legend.position='bottom')+
    scale_color_distiller(
        type = "seq",
        palette = rev('Greys'),
        direction=1)
ggrastr::rasterise(p1,dpi=500,layer='point')
dev.off()

####we'll do a reasonably low resolution here--don't want a million clusters
message("running buildSNNGraph - ", Sys.time())
snn.gr <- buildSNNGraph(sce, k = 100, use.dimred = "MNN",type='jaccard')

set.seed(100)
clust_50 <- igraph::cluster_louvain(snn.gr,resolution=1)$membership
table(clust_50)
sce$cluster<-factor(clust_50)

sce$neuron<-ifelse(sce$cluster %in% c(1,2,8,15,20,26,28),F,T)

##make boxplot by SYT1, x=cluster, y=SYT1, color=neuron v non-neuron (green v navy), no facets
pdf(file=here::here('plots','figures',
                    'supp_figures','boxplot_SYT1_cluster_coloredBySort.pdf'),h=1.5,w=4)
ggcells(sce, mapping=aes(x=cluster, y=SYT1,fill=neuron)) +
    geom_boxplot(outlier.size = 0.5)+theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 10,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_fill_manual(values=c('blue','green'))
dev.off()


##make boxplot by SYT1, x=cluster, y=detected, color=neuron v non-neuron (green v navy), no facets
pdf(file=here::here('plots','figures',
                    'supp_figures','boxplot_detectedGenes_cluster_coloredBySort.pdf'),h=1.5,w=4)
ggcells(sce, mapping=aes(x=cluster, y=detected,fill=neuron)) +
    geom_boxplot(outlier.size = 0.5)+theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 10,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_fill_manual(values=c('blue','green'))
dev.off()
###UMAPS
##cluster first
pdf(file=here::here('plots','figures','supp_figures','UMAP_initialClusters_preDropSample17.pdf'),h=6,w=6)
p<-plotUMAP(sce,text_by='cluster',colour_by='cluster',point_size=0.75,text_size=7)+theme(text=element_text(size=16))
rasterize(p,dpi=500)
dev.off()


##neurons
pdf(file=here::here('plots','figures','supp_figures','UMAP_neuronClusters_preDropSample17.pdf'),h=6,w=6)
p<-plotUMAP(sce,text_by='cluster',colour_by='neuron',point_size=0.75,text_size=7)+
    theme(text=element_text(size=16))+
    scale_color_manual(values=c('blue','#006000'))+labs(color='neuron')
rasterize(p,dpi=500)
dev.off()


##detected
pdf(file=here::here('plots','figures','supp_figures','UMAP_detectedGenes_preDropSample17.pdf'),h=6,w=6)
p<-plotUMAP(sce,text_by='cluster',colour_by='detected',point_size=0.75,text_size=7,point_alpha=1)+
    theme(text=element_text(size=16))+scale_color_distiller(
        type = "seq",
        palette = rev('Greys'),
        direction=1)+labs(color='detected genes')
rasterize(p,dpi=500)
dev.off()

###Make violin with discard
##make violin showing SYT1 vs #det genes
pdf(file=here::here('plots','figures','supp_figures','discard_round2_violin_snRNAseq.pdf'),w=4,h=3)
p1<-plotColData(sce, x = "Sample", y = "detected",
                colour_by = "discard",point_alpha=1,point_size=0.001) +
    scale_y_log10() + ylab('detected genes')+xlab('sample')+
    facet_wrap(~sce$sort,scales='free_x')+
    ggtitle("Detected genes vs round 2 discard")+
    theme(axis.text.x = element_text(size=8, angle = 90),
          legend.position='bottom')+labs(color='discard_round2')
ggrastr::rasterise(p1,dpi=500,layer='point')
dev.off()


###barplot for discard by cluster
df <- as.data.frame(colData(sce))

# Calculate proportions
df <- df %>%
    group_by(cluster, discard) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup() %>%
    group_by(cluster) %>%
    mutate(total = sum(proportion)) %>%
    mutate(proportion = proportion / total)

# Create barplot
pdf(file=here::here('plots','figures','supp_figures','barplot_discardRound2_stack.pdf'),h=6,w=8)
ggplot(df, aes(x = cluster, y = count, fill = discard)) +
    geom_bar(stat = "identity", position = "stack") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text=element_text(colour='black',size=16),
          axis.text.x = element_text(angle = 90),
          legend.position='bottom') +
    scale_fill_manual(values=c('#FF9E4A','#729ECE'))+
    ylab("nuclei") +
    xlab("cluster") +
    labs(fill = "discard round 2")
dev.off()




