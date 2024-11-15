library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(scater)
set.seed(123)

load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))
load(file=here::here('plots','spatial_palette_final.rda'))

#boxplots
#oligo
ggcells(spe, aes(x=domain, y=nmf44, fill=domain))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=spatial.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="nmf44", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        legend.position="none")
#astro
ggcells(spe, aes(x=domain, y=nmf81, fill=domain))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=spatial.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="nmf81", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        legend.position="none")
#exc postsynaptic
ggcells(spe, aes(x=domain, y=nmf13, fill=domain))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=spatial.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="nmf13", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        legend.position="none")
#inhb postsynaptic
ggcells(spe, aes(x=domain, y=nmf7, fill=domain))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=spatial.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="nmf7", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        legend.position="none")

#extended data boxplots
load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))
load(file=here::here('plots','snRNAseq_palettes.rda'))

ggcells(sce, aes(x=fine.cell.class, y=nmf13, fill=fine.cell.class))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=sn.fine.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="nmf13", y="nuclei weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        legend.position="none")
ggcells(sce, aes(x=fine.cell.class, y=nmf7, fill=fine.cell.class))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=sn.fine.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="nmf7", y="nuclei weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        legend.position="none")

#nmf13 and 7 overlap
t1 = table(cbind.data.frame(nmf13=spe$nmf13>0, nmf7=spe$nmf7>0))
t1/rowSums(t1)

ggcells(spe, aes(x=nmf13, y=nmf7, color=domain))+
  geom_point(size=.1)+scale_color_manual(values=spatial.palette)+
  facet_wrap(vars(domain))+
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none")


t2 = table(cbind.data.frame(nmf13=sce$nmf13>0, nmf7=sce$nmf7>0))
t2/rowSums(t2)

ggcells(sce, aes(x=nmf13, y=nmf7, color=fine.cell.class))+
  geom_point(size=.1)+scale_color_manual(values=sn.fine.palette)+
  scale_x_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  facet_wrap(vars(fine.cell.class))+
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none")
