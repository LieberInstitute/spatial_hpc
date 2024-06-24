###FIG7 plots
library(SingleCellExperiment)
library(scater)
library(scran)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RcppML)
library(ggspavis)
library(SpatialExperiment)
library(simplifyEnrichment)
library(pheatmap)
library(RColorBrewer)


##load sce.subset
load(file=here::here('snRNAseq_hpc','processed-data',
                     'sce','sce_subset_pats.rda'))


##load spe object
load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))

##load sce object
load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_nmf_final.rda'))

###get multicolor plotVisium()) script
source(file=here::here('code','NMF','plotvisiumRGB.R'))

##load nmf patterns
load(file=here::here('snRNAseq_hpc','processed-data','NMF','nmf_final.rda'))

###load go analysis results
load(file=here::here('snRNAseq_hpc','processed-data','NMF','go_analysis.rda'))

###get rewritten plotVisium()) script
source(file=here::here('code','NMF','plotVisium_rewrite.R'))

##########boxplots#########
sce.subset$cell.type.mouse<-sce.subset$cellType
levels(sce.subset$cell.type.mouse)[1]<-'GC'
levels(sce.subset$cell.type.mouse)[2]<-'GC'

##mouse
pdf(file=here::here('plots','figures','figure_7','mouse_nmf_10_boxplots.pdf'),h=5,w=11)
ggcells(sce.subset, mapping=aes(x=cell.type.mouse, y=nmf10,fill=condition)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 36,colour='black'))+
    theme(legend.position='bottom',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))#+scale_fill_manual(values=sn.mid.palette)
dev.off()

pdf(file=here::here('plots','figures','figure_7','mouse_nmf_14_boxplots.pdf'),h=5,w=11)
ggcells(sce.subset, mapping=aes(x=cell.type.mouse, y=nmf14,fill=condition)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 36,colour='black'))+
    theme(legend.position='bottom',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))#+scale_fill_manual(values=sn.mid.palette)
dev.off()


pdf(file=here::here('plots','figures','figure_7','mouse_nmf_13_boxplots.pdf'),h=5,w=11)
ggcells(sce.subset, mapping=aes(x=cell.type.mouse, y=nmf13,fill=condition)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 36,colour='black'))+
    theme(legend.position='bottom',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))#+scale_fill_manual(values=sn.mid.palette)
dev.off()

pdf(file=here::here('plots','figures','figure_7','mouse_nmf_91_boxplots.pdf'),h=5,w=11)
ggcells(sce.subset, mapping=aes(x=cell.type.mouse, y=nmf91,fill=condition)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 36,colour='black'))+
    theme(legend.position='bottom',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))#+scale_fill_manual(values=sn.mid.palette)
dev.off()

pdf(file=here::here('plots','figures','figure_7','mouse_nmf_20_boxplots.pdf'),h=5,w=11)
ggcells(sce.subset, mapping=aes(x=cell.type.mouse, y=PDGFB,fill=condition)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 36,colour='black'))+
    theme(legend.position='bottom',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))#+scale_fill_manual(values=sn.mid.palette)
dev.off()

##human
sce.no<-sce[,sce$fine.cell.class=='GC']
sce.no$superfine.cell.class<-droplevels(sce.no$superfine.cell.class)
pdf(file=here::here('plots','figures','figure_7','human_nmf_10_boxplots.pdf'),h=6,w=11)
ggcells(sce.no, mapping=aes(x=superfine.cell.class, y=nmf10,fill='darkgreen')) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 28,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    scale_fill_manual(values='darkgreen')
dev.off()

pdf(file=here::here('plots','figures','figure_7','human_nmf_14_boxplots.pdf'),h=6,w=11)
ggcells(sce.no, mapping=aes(x=superfine.cell.class, y=nmf14,fill='darkgreen')) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 28,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    scale_fill_manual(values='darkgreen')
dev.off()


#######multicolor spotplots#########
speb <- spe[, (colData(spe)$sample_id == 'V11L05-336_A1')]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

pdf(file=here::here('plots','figures','figure_7','rgb_gcs.pdf'),h=3.9,w=3.9)
plotVisiumRGB(speb,green ='nmf10',blue ='nmf14',image=F,highlight='neuron_cell_body',values=c('grey40','black'))
dev.off()

pdf(file=here::here('plots','figures','figure_7','rgb_nongcs.pdf'),h=3.9,w=3.9)
plotVisiumRGB(speb,yellow ='nmf91',cyan ='nmf13',pink='nmf20',image=F,highlight='neuron_cell_body',values=c('grey40','black'))
dev.off()

IQSEC2
SHANK3
KCNN3

############heatmap###########
genesplot<-c('PROX1',
            "SMAD9", "GRP", "ECE2","CNGB1","FAT4",
            'GALNT13','ZMAT4','ACVR1C','PTCHD1','TGFA',
            'JUN','FOS','ELOVL5','EGR1','NR4A1','NR4A3','RHEB',
            'SORCS3','BDNF','ETV5','RASGRF1','VGF','MICAL2','NPTX2',
             'PDZD2','NCOR2','DLGAP3','CAMKK1','SEMA5B','KCNQ2','SHANK1')

pdf(file=here::here('plots','figures','figure_7','gene_heatmap_final.pdf'),h=7,w=5)
pheatmap(x@w[genesplot,c(10,14,91,20,13)],
         cluster_rows=F,cluster_cols=F,breaks=seq(0,0.002,0.002/100),
         color=colorRampPalette(brewer.pal(n = 7, name ="Greys"))(100))
dev.off()


