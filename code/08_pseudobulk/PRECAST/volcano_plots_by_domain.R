setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
    library('SingleCellExperiment')
    library('SpatialExperiment')
    library('here')
    library('jaffelab')
    library('scater')
    library('scran')
    library('readxl')
    library('Polychrome')
    library('cluster')
    library('limma')
    library('sessioninfo')
    library('ggplot2')
    library('ggrepel')
    library('EnhancedVolcano')
    library('ggrastr')
})


options(ggrepel.max.overlaps = Inf)
pdf(file=here::here('plots','figures','supp_figures','supp_figure_volcanoes_by_domain','volcanoes.pdf'),h=5,w=5.75)
makevolcano<-function(domain){
    print(paste('printing volcano for',domain,sep=" "))
    fc<-paste0('logFC_',domain)
    fd<-paste0('fdr_',domain)
    p<-EnhancedVolcano(stats$enrichment,
                       lab = stats$enrichment$gene,
                       selectLab = stats$enrichment$gene[order(abs(stats$enrichment[[fc]]),decreasing = T)][1:8],
                       x = fc,
                       y = fd,
                       FCcutoff = 1,
                       pCutoff = 0.01,
                       ylab = "-log10 FDR",
                       legendLabels = c('n.s.','log2FC > 1','n.s',
                                        'DEGs'),
                       title = domain,
                       subtitle = NULL,
                       caption=NULL,
                       xlab='log2 fold change',
                       drawConnectors = T,
                       pointSize=0.5,
                       col=c('grey30','grey30','grey30','red2')
    )+theme( # Modify theme to remove grid lines, ticks, tick labels and set text color to black
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(color = "black"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position='none'
    )
    print(rasterize(p,dpi=350,layers = 'Point'))}
domains<- c('GCL','CA2.4','CA1','SUB','SUB.RHP','RHP','GABA','SL.SR',
            'ML','SR.SLM','SLM.SGZ','WM.1','WM.2','WM.3','Vascular','Choroid') 

sapply(domains,makevolcano)


dev.off()