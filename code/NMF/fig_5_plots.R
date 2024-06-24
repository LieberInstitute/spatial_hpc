#####Fig 5 plots
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(ggplot2)
library(schex)
library(clusterProfiler)
library(pheatmap)
library(RColorBrewer)
library(simplifyEnrichment)
library(ggspavis)

###load patterns
##load nmf patterns
load(file=here::here('snRNAseq_hpc','processed-data','NMF','nmf_final.rda'))

##load sce object
load(file=here::here('snRNAseq_hpc','processed-data',
                         'sce','sce_nmf_final.rda'))


##load palettes
load(here::here('plots','snRNAseq_palettes.rda'))

###get rewritten plotVisium()) script
source(file=here::here('code','NMF','plotVisium_rewrite.R'))


##########hexbins#############
sce_hex <- make_hexbin(sce, nbins = 100,
                       dimension_reduction = "UMAP", use_dims=c(1,2))

pdf(file=here::here('plots','figures','figure_5','nmf79_umap.pdf'),h=4,w=4)
plot_hexbin_meta(sce_hex, col="nmf79", action="median")+
    #theme(legend.position = 'bottom')+
    #legend.text = element_text(angle=90))+
    scale_fill_distiller(
        type = "seq",
        palette = rev('Oranges'),
        direction=1)+
    labs(x='UMAP1',y='UMAP2',title='nmf79 (astrocytes, ion/neurotransmitter reuptake)')+ theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
dev.off()

pdf(file=here::here('plots','figures','figure_5','nmf77_umap.pdf'),h=4,w=4)
plot_hexbin_meta(sce_hex, col="nmf77", action="median")+
    #theme(legend.position = 'bottom')+
    #legend.text = element_text(angle=90))+
    scale_fill_distiller(
        type = "seq",
        palette = rev('RdPu'),
        direction=1)+
    labs(x='UMAP1',y='UMAP2',title='nmf77 (oligodendrocytes)')+ theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
dev.off()

pdf(file=here::here('plots','figures','figure_5','nmf13_umap.pdf'),h=4,w=4)
plot_hexbin_meta(sce_hex, col="nmf13", action="median")+
    #theme(legend.position = 'bottom')+
    #legend.text = element_text(angle=90))+
    scale_fill_distiller(
        type = "seq",
        palette = rev('Greens'),
        direction=1)+
    labs(x='UMAP1',y='UMAP2',title='nmf13 (small GTPase signaling/potassium channel activity)')+ theme(
        axis.text = element_blank(),
        #axis.text.y = element_blank(),
        axis.ticks = element_blank())
dev.off()

#########boxplots#######

pdf(file=here::here('plots','figures','figure_5','nmf_77_boxplots.pdf'),h=2,w=3)
ggcells(sce, mapping=aes(x=fine.cell.class, y=nmf77,fill=mid.cell.class)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 9,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_fill_manual(values=sn.mid.palette)
dev.off()

pdf(file=here::here('plots','figures','figure_5','nmf_79_boxplots.pdf'),h=2,w=3)
ggcells(sce, mapping=aes(x=fine.cell.class, y=nmf79,fill=mid.cell.class)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 9,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_fill_manual(values=sn.mid.palette)
dev.off()

pdf(file=here::here('plots','figures','figure_5','nmf_13_boxplots.pdf'),h=2,w=3)
ggcells(sce, mapping=aes(x=fine.cell.class, y=nmf13,fill=mid.cell.class)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 9,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_fill_manual(values=sn.mid.palette)
dev.off()
##########heatmap###########
###marker genes
features<-c('PLP1','CNP','SCD','CRYAB','MAG',
            'SLC1A2','SLC1A3','SLC25A18','NTM','GPC5',
            'PDZD2','DLGAP3','PDE2A','KCNQ2','RASAL1')
heat<-x@w[features,c('nmf77','nmf79','nmf13')]

pdf(file=here::here('plots','figures','figure_5','marker_gene_heatmap.pdf'),h=5,w=3)
pheatmap(heat,cluster_rows=F,cluster_cols=F,
         color=colorRampPalette(brewer.pal(n = 7, name ="Greys"))(100))
dev.off()

#########ORA dotplot##########
##load go analysis
load(file=here::here('snRNAseq_hpc','processed-data','NMF','go_analysis.rda'))

##subset to patterns of interest
go<-go[c('nmf77','nmf79','nmf13')]


##get top terms
top1<-go[[1]]$Description[order(go[[1]]$p.adjust)][1:5]
index1<-which(go[[1]]$Description %in% top1)
go[[1]]<-subset_enrichResult(go[[1]],index1)

##subsetting manually to exclude redundant terms
top2<-c('acyl-CoA metabolic process','organic acid catabolic process',
        'L-aspartate transmembrane transport','amino acid metabolic process',
        'L-glutamate transmembrane transport')
index2<-which(go[[2]]$Description %in% top2)
go[[2]]<-subset_enrichResult(go[[2]],index2)

top3<-go[[3]]$Description[order(go[[3]]$p.adjust)][1:5]
index3<-which(go[[3]]$Description %in% top3)
go[[3]]<-subset_enrichResult(go[[3]],index3)

merged<-merge_result(go)

pdf(file=here::here('plots','figures','figure_5','dotplot.pdf'),h=10,w=5.4)
dotplot(merged,showCategory=5,includeAll=F)+scale_fill_distiller(
    type = "seq",
    palette = rev('Greys'),
    direction=-1)
dev.off()

############spotplots############
##load data
load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))

##set up spotplots
brains <- unique(spe$brnum)


speb <- spe[, (colData(spe)$brnum == brains[[9]])]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)


##make background palette
spatial.palette3<-c("#BEDDBA", "#eae8e4", "#E8BBC6","#A1BAD8")
names(spatial.palette3)<-levels(speb$broad.domain)

pdf(file=here::here('plots','figures','figure_5','nmf_77_spotplots.pdf'),h=8,w=8)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'nmf77',
    highlight = "broad.domain",
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    #palette = spatial.palette,
    image=F,
    values=spatial.palette3)+
    scale_fill_distiller(
        type = "seq",
        palette = rev('Greys'),
        direction=1)
dev.off()

pdf(file=here::here('plots','figures','figure_5','nmf_79_spotplots.pdf'),h=8,w=8)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'nmf79',
    highlight = "broad.domain",
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    #palette = spatial.palette,
    image=F,
    values=spatial.palette3)+
    scale_fill_distiller(
        type = "seq",
        palette = rev('Greys'),
        direction=1)
dev.off()

pdf(file=here::here('plots','figures','figure_5','nmf_13_spotplots.pdf'),h=8,w=8)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'nmf13',
    highlight = "broad.domain",
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    #palette = spatial.palette,
    image=F,
    values=spatial.palette3)+
    scale_fill_distiller(
        type = "seq",
        palette = rev('Greys'),
        direction=1)
dev.off()
