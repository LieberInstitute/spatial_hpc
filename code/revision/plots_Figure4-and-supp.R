library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(scater)
set.seed(123)

load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))
load(file=here::here('plots','spatial_palette_final.rda'))

#oligo
p1 <- ggcells(spe, aes(x=domain, y=nmf44, fill=domain))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=spatial.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Oligodendrocytes (nmf44)", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())
#astro
p2 <- ggcells(spe, aes(x=domain, y=nmf81, fill=domain))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=spatial.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Astrocytes (nmf81)", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())
#exc postsynaptic
p3 <- ggcells(spe, aes(x=domain, y=nmf13, fill=domain))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=spatial.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Excitatory postsynaptic (nmf13)", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())
#inhb postsynaptic
p4 <- ggcells(spe, aes(x=domain, y=nmf7, fill=domain))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=spatial.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Inhibitory postsynaptic (nmf7)", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())


#extended data boxplots
load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))
load(file=here::here('plots','snRNAseq_palettes.rda'))

#PrsPas --> L2/3.1; PrsENT --> L2/3.2
sce$fine.cell.class <- factor(sce$fine.cell.class, levels=levels(sce$fine.cell.class),
                              labels=c(levels(sce$fine.cell.class)[1:7],
                                       "L2/3.1","L2/3.2",
                                       levels(sce$fine.cell.class)[10:24]))
names(sn.fine.palette) <- levels(sce$fine.cell.class)

p5 <- ggcells(sce, aes(x=fine.cell.class, y=nmf13, fill=fine.cell.class))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=sn.fine.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Excitatory postsynaptic (nmf13) - snRNAseq", y="nuclei weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=6),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())
p6 <- ggcells(sce, aes(x=fine.cell.class, y=nmf7, fill=fine.cell.class))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=sn.fine.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Inhibitory postsynaptic (nmf7) - snRNAseq", y="nuclei weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=6),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())

#nmf13 and 7 overlap

p7 <- ggcells(spe, aes(x=nmf13, y=nmf7, color=domain))+
  geom_point(size=.1)+scale_color_manual(values=spatial.palette)+
  scale_x_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  facet_wrap(vars(domain))+ggtitle("Co-occurrence of nmf13 and nmf7 in SRT")+
  theme(text=element_text(size=10), 
        axis.text.x=element_text(angle=45, hjust=1, size=6), 
        axis.text.y=element_text(size=6), legend.position="none")

p8 <- ggcells(sce, aes(x=nmf13, y=nmf7, color=fine.cell.class))+
  geom_point(size=.1)+scale_color_manual(values=sn.fine.palette)+
  scale_x_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  facet_wrap(vars(fine.cell.class))+ggtitle("Co-occurrence of nmf13 and nmf7 in snRNAseq")+
  theme(text=element_text(size=10), panel.spacing=unit(3,'points'), 
        strip.text=element_text(margin=margin(2,2,2,2, unit="pt")),
        axis.text.x=element_text(angle=45, hjust=1, size=6),
        axis.text.y=element_text(size=6), legend.position="none")

laymat = rbind(c(1,2),c(1,2),c(3,4),c(3,4),c(5,6),c(5,6),c(7,8),c(7,8),c(7,8),c(7,8))

gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, layout_matrix=laymat)

pdf(file = "plots/revision/Figure4_supp-boxplots-scatter.pdf",
    width=8.5, height=11)
gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, layout_matrix=laymat)
dev.off()

