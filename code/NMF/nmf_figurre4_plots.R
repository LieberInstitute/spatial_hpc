library(scales)
library(SpatialExperiment)
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggspavis)

#make mossy_cell logical plot
pdf('nmf52.pdf',h=5,w=6)
plot<-plotUMAP(sce,colour_by='nmf52',point_size=0.01)+
    scale_color_distiller(
         type = "seq",
         palette = rev('RdPu'),
         direction=1) +
            labs(color='nmf52')
ggrastr::rasterize(plot,layer='point',dpi=500)
dev.off()
pdf('mossycells.pdf',h=5,w=6)
plot<-plotUMAP(sce,colour_by='mossy_cell',point_size=0.01)+
    scale_color_manual(values=c('#feebe2','#ae017e'))+
    labs(color='mossy cells')
ggrastr::rasterize(plot,layer='point',dpi=500)
dev.off()

####Figure 5
###Heatmap, snRNA-seq clusters
###Heatmap, visium clusters
###Violins,
###Spotplot #1
pdf(file=here::here('plots/figures/figure_5/spotplot1.pdf'),h=7,w=7)
plotVisiumRGB(speb,vars=c('nmf11','nmf15','nmf52','nmf61'),image=F,highlight='neuron')
dev.off()
##Spotplot #2
pdf(file=here::here('plots/figures/figure_5/spotplot2.pdf'),h=7,w=7)
plotVisiumRGB(speb,vars=c('nmf40','nmf17','nmf32','nmf54'),image=F,highlight='neuron')
dev.off()
# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="https://dec2021.archive.ensembl.org/")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="https://dec2021.archive.ensembl.org/")
    genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    humanx <- genesV2
    # Print the first 6 genes found to the screen
    print(head(humanx))
    return(humanx)
}

genes<-genes[genes$MGI.symbol %in% rownames(results@featureLoadings),]
genes<-genes[genes$HGNC.symbol %in% rownames(sce),]
genes<-genes[!duplicated(genes$MGI.symbol),]
genes<-genes[!duplicated(genes$HGNC.symbol),]


i<-intersect(genes$HGNC.symbol,rownames(sce))
j<-intersect(genes$MGI.symbol,rownames(results@featureLoadings))

#loadings<-x_mouse@w
loadings<-results@featureLoadings
loadings<-loadings[rownames(loadings) %in% j,]
loadings<-loadings[match(genes$MGI.symbol,rownames(loadings)),]

rownames(loadings)<-genes$HGNC.symbol


sce2<-sce[rownames(sce) %in% i,]
data<-as.matrix(logcounts(sce2))

proj<-projectR(
    data=data,
    loadings=loadings,
    full = FALSE,
    family = "gaussianff",
    bootstrapPval = FALSE,
    bootIter = 1000
)

##NMF figures:
##Fig 4B
pdf(file=here::here('plots','figures','figure_4','pat5_umap.pdf'),h=4,w=4)
plot<-plotUMAP(sce,colour_by='nmf5', point_size=0.25)+scale_color_distiller(
    type = "seq",
    palette = rev('PuRd'),
    direction=1)+labs(color='Pattern 15')#+theme(legend.position='bottom')
ggrastr::rasterize(plot,layer='point',dpi=500)
dev.off()

pdf(file=here::here('plots','figures','figure_4','pat81_umap.pdf'),h=4,w=4)
plot<-plotUMAP(sce,colour_by='nmf81', point_size=0.25)+scale_color_distiller(
    type = "seq",
    palette = rev('PuRd'),
    direction=1)+labs(color='Pattern 81 (astro-specific)')+theme(legend.position='bottom')
ggrastr::rasterize(plot,layer='point',dpi=500)
dev.off()


pdf(file=here::here('plots','figures','figure_4','pat6_umap.pdf'),h=4,w=4)
plot<-plotUMAP(sce,colour_by='nmf6', point_size=0.25)+scale_color_distiller(
    type = "seq",
    palette = rev('PuRd'),
    direction=1)+labs(color='Pattern 6 (neuron-specific)')+theme(legend.position='bottom')
ggrastr::rasterize(plot,layer='point',dpi=500)
dev.off()

##Figure 4c
pdf(file=here::here('plots','figures','figure_4','pat5_boxplot.pdf'),h=4,w=4)
ggcells(sce, mapping=aes(x=cell.class, y=nmf5,fill=broad.class)) +
    geom_boxplot()+theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
dev.off()

pdf(file=here::here('plots','figures','figure_4','pat81_boxplot.pdf'),h=4,w=4)
ggcells(sce, mapping=aes(x=cell.class, y=nmf81,fill=broad.class)) +
    geom_boxplot()+theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
dev.off()

##Figure 4d
speb <- spe[, (colData(spe)$sample_id == 'V11L05-333_D1')]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

##
plotVisium(
    speb,
    spots = TRUE,
    fill = 'COX6A1',
    #    highlight = 'neuron',
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    #palette = as.vector(palette36.colors(17)),
    image=F
)+scale_fill_distiller(
    type = "seq",
    palette = rev('RdPu'),
    direction=1)



