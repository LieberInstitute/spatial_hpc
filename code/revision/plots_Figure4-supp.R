library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(scater)
set.seed(123)

load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))
load(file=here::here('plots','spatial_palette_final.rda'))

##### domain

#oligo
p1 <- ggcells(spe, aes(x=domain, y=nmf44, fill=domain))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=spatial.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Oligodendrocytes (nmf44)", subtitle="Spatial domain", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())
#astro
p2 <- ggcells(spe, aes(x=domain, y=nmf81, fill=domain))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=spatial.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Astrocytes (nmf81)", subtitle="Spatial domain", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())

#check against astro 79 and oligo 77
p1.1 <- ggcells(spe, aes(x=domain, y=nmf77, fill=domain))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=spatial.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Oligodendrocytes (nmf77)", subtitle="Spatial domain", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())

p1.2 <- ggcells(spe, aes(x=domain, y=nmf38, fill=domain))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=spatial.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Oligodendrocytes (nmf38)", subtitle="Spatial domain", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())

p2.1 <- ggcells(spe, aes(x=domain, y=nmf79, fill=domain))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=spatial.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Astrocytes (nmf79)", subtitle="Spatial domain", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())

##### mannual annotation

load("processed-data/manual_annotation_csv/compiled_annotation_before_match.Rdata")
test = read.csv("processed-data/manual_annotation_csv/spatialLIBD_ManualAnnotation_2023-04-12_Br2720_all.csv") %>%
  mutate(sample_id= as.character(factor(sample_id, levels=c("Br2720_A1","Br2720_B1","Br2720_C1","Br2720_D1"),
                                        labels=c("V12F14-051_A1","V12F14-051_B1","V12F14-051_C1","V12F14-051_D1"))),
         spot_name=paste(spot_name, sample_id, sep="_"))
csv2 = rbind(csv2, test)

dim(spe) #31483 150917
spe$spot_id = paste(sapply(strsplit(rownames(colData(spe)), "_"), function(x) x[[1]]), spe$slide, spe$array, sep="_")
spe.manual = spe[,spe$spot_id %in% csv2$spot_name]
dim(spe.manual) #31483 150520

rownames(csv2) = csv2$spot_name
csv2 = csv2[spe.manual$spot_id,]
identical(spe.manual$spot_id, csv2$spot_name)

spe.manual$ManualAnnotation = factor(csv2$ManualAnnotation, 
                                     levels=c("GCL","CA4","PCL-CA3","PCL-CA1","SUB","CTX","THAL",
                                              "SL","SR","ML","SGZ","SLM","SO","WM","CP"))


man.pal = c("THAL"="#1e1eff","CTX"="#5ffffb", "SUB"="#add294", "PCL-CA1"="#00dc00",
            "PCL-CA3"="#00a000", "CA4"="#B0BF1A",
            "GCL"="#005000", "SGZ"="#dfa56e", "ML"="#c1c1c1", "SL"="#444444", "SR"="#828E84", "SLM"="tan4",
            "SO"="#A698AE", "WM"="#ff3ffc", "CP"="#00006a")

#oligo
p3 <- ggcells(spe.manual, aes(x=ManualAnnotation, y=nmf44, fill=ManualAnnotation))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=man.pal)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Oligodendrocytes (nmf44)", subtitle="Manual annotation", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())
#astro
p4 <- ggcells(spe.manual, aes(x=ManualAnnotation, y=nmf81, fill=ManualAnnotation))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=man.pal)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Astrocytes (nmf81)", subtitle="Manual annotation", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())


##### snRNAseq
load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))
load(file=here::here('plots','snRNAseq_palettes.rda'))

sce$fine.cell.class = factor(sce$fine.cell.class, levels=levels(sce$fine.cell.class),
                             labels=c(levels(sce$fine.cell.class)[1:7], "L2/3.1","L2/3.2",
                                      levels(sce$fine.cell.class)[10:24]))

names(sn.fine.palette) = levels(sce$fine.cell.class)

#olig
p5 <- ggcells(sce, aes(x=fine.cell.class, y=nmf44, fill=fine.cell.class))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=sn.fine.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Oligodendrocytes (nmf44)", subtitle="snRNA-seq", y="nuclei weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())
#astro
p6 <- ggcells(sce, aes(x=fine.cell.class, y=nmf81, fill=fine.cell.class))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=sn.fine.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Astrocytes (nmf81)", subtitle="snRNA-seq", y="nuclei weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())


#check against oligo 77
p5.1 <- ggcells(sce, aes(x=fine.cell.class, y=nmf77, fill=fine.cell.class))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=sn.fine.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Oligodendrocytes (nmf77)", subtitle="snRNA-seq", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())


##### save astro and oligo

pdf(file = "plots/revision/Figure4_astro-oligo-boxplots.pdf",
    width=8.5, height=9)
gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, ncol=2)
dev.off()




# now do a different supp figure for trans-neuronal




#exc postsynaptic
p1 <- ggcells(spe, aes(x=domain, y=nmf13, fill=domain))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=spatial.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Excitatory postsynaptic (nmf13)", subtitle="Spatial domain", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())
#inhb postsynaptic
p2 <- ggcells(spe, aes(x=domain, y=nmf7, fill=domain))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=spatial.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Inhibitory postsynaptic (nmf7)", subtitle="Spatial domain", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())
#exc postsynaptic
p3 <- ggcells(spe.manual, aes(x=ManualAnnotation, y=nmf13, fill=ManualAnnotation))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=man.pal)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Excitatory postsynaptic (nmf13)", subtitle="Manual annotation", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())
#inhb postsynaptic
p4 <- ggcells(spe.manual, aes(x=ManualAnnotation, y=nmf7, fill=ManualAnnotation))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=man.pal)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Inhibitory postsynaptic (nmf7)", subtitle="Manual annotation", y="spot weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())
#exc postsynaptic
p5 <- ggcells(sce, aes(x=fine.cell.class, y=nmf13, fill=fine.cell.class))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=sn.fine.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Excitatory postsynaptic (nmf13)", subtitle="snRNA-seq", y="nuclei weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=6),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())
#inhb postsynaptic
p6 <- ggcells(sce, aes(x=fine.cell.class, y=nmf7, fill=fine.cell.class))+
  geom_boxplot(outlier.size=.1)+theme_bw()+
  scale_fill_manual(values=sn.fine.palette)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  labs(title="Inhibitory postsynaptic (nmf7)", subtitle="snRNA-seq", y="nuclei weights")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=6),
        axis.text.y=element_text(size=8),
        legend.position="none", axis.title.x=element_blank())

#nmf13 and 7 overlap

p7 <- ggcells(spe, aes(x=nmf13, y=nmf7, color=domain))+
  geom_point(size=.1)+scale_color_manual(values=spatial.palette)+
  scale_x_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  facet_wrap(vars(domain))+
  labs(title="Co-occurrence of nmf13 and nmf7", subtitle="Spatial domain")+
  theme(text=element_text(size=10), 
        axis.text.x=element_text(angle=45, hjust=1, size=6), 
        axis.text.y=element_text(size=6), legend.position="none")

p8 <- ggcells(spe.manual, aes(x=nmf13, y=nmf7, color=ManualAnnotation))+
  geom_point(size=.1)+scale_color_manual(values=man.pal)+
  scale_x_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  facet_wrap(vars(ManualAnnotation))+
  labs(title="Co-occurrence of nmf13 and nmf7", subtitle="Manual annotation")+
  theme(text=element_text(size=10), 
        axis.text.x=element_text(angle=45, hjust=1, size=6), 
        axis.text.y=element_text(size=6), legend.position="none")

p9 <- ggcells(sce, aes(x=nmf13, y=nmf7, color=fine.cell.class))+
  geom_point(size=.1)+scale_color_manual(values=sn.fine.palette)+
  scale_x_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
  facet_wrap(vars(fine.cell.class))+
  labs(title="Co-occurrence of nmf13 and nmf7", subtitle="snRNA-seq")+
  theme(text=element_text(size=10), panel.spacing=unit(3,'points'), 
        strip.text=element_text(margin=margin(2,2,2,2, unit="pt")),
        axis.text.x=element_text(angle=45, hjust=1, size=6),
        axis.text.y=element_text(size=6), legend.position="none")


laymat = rbind(c(1,2),c(1,2),c(3,4),c(3,4),c(5,6),c(5,6),c(7,8),c(7,8),c(7,8),c(7,8))


pdf(file = "plots/revision/Figure4_trans-nrn-boxplots-scatter.pdf",
    width=8.5, height=11.5)
gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, p7, p9, layout_matrix=laymat)
dev.off()
