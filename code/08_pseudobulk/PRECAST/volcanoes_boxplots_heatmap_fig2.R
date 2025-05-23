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
})

## load DE data
load(file=here::here('processed-data','08_pseudobulk','PRECAST','visiumHE_DE_stats_tissueLevel.rda'))
load(file=here::here('processed-data','08_pseudobulk','PRECAST','spe_pseudo_HE.rda'))



thresh_fdr <- 0.01
thresh_logfc <- log2(2)
fdrs_gene_ids <- rowData(spe_pseudo)$gene_id
fdrs_gene_names <- rowData(spe_pseudo)$gene_name

pdf(file=here::here('plots','figures','figure_2','neuropil_volcano.pdf'),h=5,w=5.75)
    p<-EnhancedVolcano(stats$enrichment,
                    lab = stats$enrichment$gene,
                    selectLab = c('MT-ND4','ETNPPL','SLC1A3','MT-ND5','SULF1','ATP1B2'),
                    x = 'logFC_Neuropil',
                    y = 'fdr_Neuropil',
                    FCcutoff = 1,
                    pCutoff = 0.01,
                    ylab = "-log10 FDR",
                    legendLabels = c('n.s.','log2FC > 1','n.s',
                                     'DEGs'),
                    title = 'Neuropil',
                    subtitle = NULL,
                    caption=NULL,
                    xlab='log2 fold change',
                    drawConnectors = T,
                    pointSize=1,
                    col=c('grey30','grey30','grey30','red2')
    )+theme( # Modify theme to remove grid lines, ticks, tick labels and set text color to black
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(color = "black"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position='none'
    )
rasterize(p,dpi=350,layers = 'Point')

dev.off()

pdf(file=here::here('plots','figures','figure_2','neuron_volcano.pdf'),h=5,w=5.75)
p<-EnhancedVolcano(stats$enrichment,
                lab = stats$enrichment$gene,
                selectLab = c('PREPL','CLSTN3','ENO2',
                              'RPL13A','RPS6','IFITM3'),
                max.overlaps=100,
                x = 'logFC_Neuron',
                y = 'fdr_Neuron',
                FCcutoff = 1,
                pCutoff = 0.01,
                ylab = "-log10 FDR",
                legendLabels = c('n.s.','n.s.','n.s.',
                                 'DEGs'),
                col=c('grey30','grey30','grey30','red2'),
                title = 'Neuron',
                subtitle = NULL,
                caption=NULL,
                xlab='log2 fold change',
                drawConnectors = T,
                pointSize=1
)+theme( # Modify theme to remove grid lines, ticks, tick labels and set text color to black
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = 'none'
)
rasterize(p,dpi=350,layers = 'Point')
dev.off()

pdf(file=here::here('plots','figures','figure_2','wm_volcano.pdf'),h=5,w=5.75)
p<-EnhancedVolcano(stats$enrichment,
                lab = stats$enrichment$gene,
                selectLab = c('ABCA2','SHTN1','TMEM165','DBNDD2',
                              'PTMS','DKK3'),
                max.overlaps=100,
                x = 'logFC_WM',
                y = 'fdr_WM',
                FCcutoff = 1,
                pCutoff = 0.01,
                ylab = "-log10 FDR",
                legendLabels = c('n.s.','n.s.','n.s.',
                                 'DEGs'),
                col=c('grey30','grey30','grey30','red2'),
                title = 'WM',
                subtitle = NULL,
                caption=NULL,
                xlab='log2 fold change',
                drawConnectors = T,
                pointSize=1
)+theme( # Modify theme to remove grid lines, ticks, tick labels and set text color to black
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = 'none'
)
rasterize(p,dpi=350,layers = 'Point')
dev.off()


pdf(file=here::here('plots','figures','figure_2','vasc_volcano.pdf'),h=5,w=5.75)
p<-EnhancedVolcano(stats$enrichment,
                lab = stats$enrichment$gene,
                selectLab = c('TPM2','TAGLN','FLNA','MYH11',
                              'KIF5C','ANK3'),
                max.overlaps=100,
                x = 'logFC_Vasc_CSF',
                y = 'fdr_Vasc_CSF',
                FCcutoff = 1,
                pCutoff = 0.01,
                ylab = "-log10 FDR",
                legendLabels = c('n.s.','n.s.','n.s.',
                                 'DEGs'),
                col=c('grey30','grey30','grey30','red2'),
                title = 'Vasc_CSF',
                subtitle = NULL,
                caption=NULL,
                xlab='log2 fold change',
                drawConnectors = T,
                pointSize=1
)+theme( # Modify theme to remove grid lines, ticks, tick labels and set text color to black
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position='none'
)
rasterize(p,dpi=350,layers = 'Point')
dev.off()


broad_cols<-spatial.palette2
pdf(file=here::here('plots','figures','figure_2','clstn3.pdf'),w=4.2,h=4)
ggcells(spe_pseudo, mapping=aes(x=domain, y=CLSTN3,fill=broad.domain)) +
    geom_boxplot()+theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=broad_cols)
dev.off()

pdf(file=here::here('plots','figures','figure_2','SLC1A3.pdf'),w=4.2,h=4)
ggcells(spe_pseudo, mapping=aes(x=domain, y=SLC1A3,fill=broad.domain)) +
    geom_boxplot()+theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    scale_fill_manual(values=broad_cols)
dev.off()

pdf(file=here::here('plots','figures','figure_2','SHTN1.pdf'),w=4.2,h=4)
ggcells(spe_pseudo, mapping=aes(x=domain, y=SHTN1,fill=broad.domain)) +
    geom_boxplot()+theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    scale_fill_manual(values=broad_cols)
dev.off()

pdf(file=here::here('plots','figures','figure_2','TPM2.pdf'),w=4.2,h=4)
ggcells(spe_pseudo, mapping=aes(x=domain, y=TPM2,fill=broad.domain)) +
    geom_boxplot()+theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    scale_fill_manual(values=broad_cols)
dev.off()

##Heatmap
features <- c("PPFIA2", "AMPH", "KCNG2", "PRKCG", "FNDC1",
              "GFRA1", "TOX", "SLC17A6", "MEF2C", "GAD2", "MIF",
              "FABP7", "APOC1", "MT1G", "MAN1A2", "NTRK2", "SFRP2",
              "ABCA2", "MOBP", "MTURN", "PHLDB1", "MFAP4", "ACTA2", "PRLR")
pdf(file=here::here('plots','figures','figure_2','cluster_markers_heatmap.pdf'),h=3.7,w=4)
    plotGroupedHeatmap(spe_pseudo,group = 'domain',features=features,
                   scale=T,center=T,cluster_rows=F,cluster_cols=F,
                   exprs_values = 'logcounts',color=viridis::magma(25))
dev.off()
