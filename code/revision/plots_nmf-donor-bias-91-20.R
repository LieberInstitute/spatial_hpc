library(SingleCellExperiment)
library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scater)

set.seed(123)

load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))
avg.sce = rowMeans(logcounts(sce))

load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))

#load nmfs
load(file=here::here('snRNAseq_hpc','processed-data','NMF','nmf_final.rda'))
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]
avg.sce = avg.sce[rownames(loads)]

#distribution across broad domains
table(filter(as.data.frame(colData(spe)), nmf91>0)$broad.domain)
(4011+2451)/(4011+2451+1847+1071)
(4011)/(4011+2451+1847+1071)
table(filter(as.data.frame(colData(spe)), nmf20>0)$broad.domain)
(1919+147)/(1919+147+26+17)
1919/(1919+147+26+17)

#pie
p1 <- ggplot(filter(as.data.frame(colData(sce)), nmf91>0),
       aes(x="", fill=brnum))+
  geom_bar(stat="count", position="fill", color="black")+
  scale_fill_manual(values=c("#313695",#2720
                             "#33a02c",#2743
                             "#fb9a99",#3492
                             "#ff7f00",#6423
                             "#6a3d9a",#6432
                             "grey",#6471
                             "#ffff99",#6522
                             "#8dd3c7",#8325
                             "#b2182b",#8492
                             "grey40" #8667
  ))+labs(fill="Donor", title="snRNAseq: nmf91")+
  coord_polar(theta="y")+theme_void()+
  theme(plot.title=element_text(size=18), aspect.ratio=1)

p2 <- ggplot(filter(as.data.frame(colData(spe)), nmf91>0),
       aes(x="", fill= factor(brnum, levels=levels(spe$brnum)[c(10,1:9)])))+
  geom_bar(stat="count", position="fill", color="black")+
  scale_fill_manual(values=c("#313695",#2720
                             "#33a02c",#2743
                             "#fb9a99",#3492
                             "#ff7f00",#6423
                             "#6a3d9a",#6432
                             "grey",#6471
                             "#ffff99",#6522
                             "#8dd3c7",#8325
                             "#b2182b",#8492
                             "grey40" #8667
  ))+labs(fill="Donor", title="SRT: nmf91")+
  coord_polar(theta="y")+theme_void()+
  theme(plot.title=element_text(size=18), aspect.ratio=1)

p3 <- ggplot(filter(as.data.frame(colData(sce)), nmf20>0),
       aes(x="", fill=brnum))+
  geom_bar(stat="count", position="fill", color="black")+
  scale_fill_manual(values=c("#313695",#2720
                             "#33a02c",#2743
                             "#fb9a99",#3492
                             "#ff7f00",#6423
                             "#6a3d9a",#6432
                             "grey",#6471
                             "#ffff99",#6522
                             "#8dd3c7",#8325
                             "#b2182b",#8492
                             "grey40" #8667
  ))+labs(fill="Donor", title="snRNAseq: nmf20")+
  coord_polar(theta="y")+theme_void()+
  theme(plot.title=element_text(size=18), aspect.ratio=1)

p4 <- ggplot(filter(as.data.frame(colData(spe)), nmf20>0),
       aes(x="", fill= factor(brnum, levels=levels(spe$brnum)[c(10,1:9)])))+
  geom_bar(stat="count", position="fill", color="black")+
  scale_fill_manual(values=c("#313695",#2720
                             "#33a02c",#2743
                             "#fb9a99",#3492
                             "#ff7f00",#6423
                             "#6a3d9a",#6432
                             "grey",#6471
                             "#ffff99",#6522
                             "#8dd3c7",#8325
                             "#b2182b",#8492
                             "grey40" #8667
  ))+labs(fill="Donor", title="SRT: nmf20")+
  coord_polar(theta="y")+theme_void()+
  theme(plot.title=element_text(size=18), aspect.ratio=1)

ggsave(file="plots/revision/ED_nmf91-nmf20_pie-charts.pdf", 
       gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2),
       height=8, width=8, units="in")

#weight ~ avg expr scatter
nmf91.df = tibble::rownames_to_column(cbind.data.frame(nmf91=loads[,"nmf91"], avg.expr=avg.sce), var="gene") %>%
  mutate(is_ieg=gene %in% c("NR4A3","JUND","JUN","FOS","NR4A1","NR4A2","EGR1","JUNB","FOSL2"))
p5 <- ggplot(nmf91.df, aes(x=avg.expr, y=nmf91, color=is_ieg))+
  geom_point()+
  geom_point(data=filter(nmf91.df, is_ieg==TRUE))+
  geom_text_repel(data=filter(nmf91.df, is_ieg==TRUE), aes(label=gene),
                  min.segment.length = 0.3, #point.padding=.5, 
                  box.padding=.5, color="black", seed=123)+
  theme_bw()+scale_color_manual(values=c("grey","red"))+
  labs(x="logcount expr. (all nuclei)", y="gene weight", title="nmf91")+
  theme(plot.title=element_text(size=16))

nmf20.sort = sort(loads[,"nmf20"], decreasing=T)
nmf20.df = tibble::rownames_to_column(cbind.data.frame(nmf20=loads[,"nmf20"], avg.expr=avg.sce), var="gene") %>%
  mutate(is_ieg=gene %in% names(nmf20.sort)[1:3])
p6 <- ggplot(nmf20.df, aes(x=avg.expr, y=nmf20, color=is_ieg))+
  geom_point()+
  geom_point(data=filter(nmf20.df, is_ieg==TRUE))+
  geom_text_repel(data=filter(nmf20.df, is_ieg==TRUE), aes(label=gene),
                  min.segment.length = 0.3, #point.padding=.5, 
                  box.padding=.5, color="black", seed=123)+
  theme_bw()+scale_color_manual(values=c("grey","red"))+
  labs(x="logcount expr. (all nuclei)", y="gene weight", title="nmf20")+
  theme(plot.title=element_text(size=16))

ggsave(file="plots/revision/ED_nmf91-nmf20_gene-weight-scatter.pdf", 
       gridExtra::grid.arrange(ggrastr::rasterize(p5,layer='point',dpi=350), ggrastr::rasterize(p6,layer='point',dpi=350),
                               ncol=2),
       height=5, width=10, units="in")

#settled on dotplots instead of boxplots
sce$brnum_nrn = ifelse(sce$broad.cell.class=="Neuron", paste(sce$brnum, "nrn", sep="_"),
                       paste(sce$brnum, "other", sep="_"))
sce$brnum_nrn = factor(sce$brnum_nrn, levels=rev(paste0(rep(levels(sce$brnum), each=2), rep(c("_nrn","_other"), length.out=20))))

dplot1 <- plotDots(sce, features=names(nmf20.sort)[3:1], group = "brnum_nrn")+
  scale_color_gradient(low="white", high="black", limits=c(0,3.1))+
  scale_size(limits=c(0,1))+coord_flip()+
  labs(size="Prop. nuclei",color="Avg. expr.")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=12), axis.title.y=element_blank(),
        text=element_text(size=10), axis.title.x=element_blank())

dplot2 <- plotDots(sce, features=c("JUN","FOS","NR4A1"), group = "brnum_nrn")+
  scale_color_gradient(low="white", high="black", limits=c(0,1.3))+
  scale_size(limits=c(0,1))+coord_flip()+
  labs(size="Prop. nuclei",color="Avg. expr.")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=12), axis.title.y=element_blank(),
        text=element_text(size=10), axis.title.x=element_blank())

spe$brnum_nrn = ifelse(spe$broad.domain=="Neuron", paste(spe$brnum, "nrn", sep="_"),
                         paste(spe$brnum, "other", sep="_"))
spe$brnum_nrn = factor(spe$brnum_nrn, levels=rev(paste0(rep(levels(sce$brnum), each=2), rep(c("_nrn","_other"), length.out=20))))
dplot3 <- plotDots(spe, features=names(nmf20.sort)[3:1], group = "brnum_nrn")+
  scale_color_gradient(low="white", high="black", limits=c(0,.3))+
  scale_size(limits=c(0,1))+coord_flip()+
  labs(size="Prop. spots",color="Avg. expr.")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y=element_text(size=12), axis.title.y=element_blank(),
        text=element_text(size=10), axis.title.x=element_blank())
dplot4 <- plotDots(spe, features=c("JUN","FOS","NR4A1"), group = "brnum_nrn")+
    scale_color_gradient(low="white", high="black", limits=c(0,1))+
    scale_size(limits=c(0,1))+coord_flip()+
    labs(size="Prop. spots",color="Avg. expr.")+
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
          axis.text.y=element_text(size=12), axis.title.y=element_blank(),
          text=element_text(size=10), axis.title.x=element_blank())

ggsave(file="plots/revision/ED_nmf91-nmf20_top-gene-dotplot.pdf",
       gridExtra::grid.arrange(dplot1, dplot3, dplot2, dplot4, ncol=4),
       height=6, width=12, units="in"
)

#boxplots
ggplot(mutate(as.data.frame(colData(sce)), is_nrn=broad.cell.class=="Neuron",
              ieg=logcounts(sce)["JUN",]),
       aes(x=brnum, y=ieg, fill=is_nrn))+
  geom_boxplot(outlier.size=.1)+
  theme_bw()+
  scale_fill_manual(values=c("grey50","red"))+
  labs(y="logcounts expr.", x="", title="JUN (snRNAseq)")

ggplot(mutate(as.data.frame(colData(spe)), is_nrn=broad.domain=="Neuron",
              ieg=logcounts(spe)["JUN",]),
       aes(x=factor(brnum, levels=levels(spe$brnum)[c(10,1:9)]), y=ieg, fill=is_nrn))+
  geom_boxplot(outlier.size=.1)+
  theme_bw()+
  scale_fill_manual(values=c("grey50","red"))+
  labs(y="logcounts expr.", x="", title="JUN (SRT)")


ggplot(mutate(as.data.frame(colData(sce)), is_nrn=broad.cell.class=="Neuron",
              srg=logcounts(sce)["HOMER1",]),
       aes(x=brnum, y=srg, fill=is_nrn))+
  geom_boxplot(outlier.size=.1)+
  theme_bw()+
  scale_fill_manual(values=c("grey50","red"))+
  labs(y="logcounts expr.", x="", title="HOMER1 (snRNAseq)")

ggplot(mutate(as.data.frame(colData(spe)), is_nrn=broad.domain=="Neuron",
              srg=logcounts(spe)["HOMER1",]),
       aes(x=factor(brnum, levels=levels(spe$brnum)[c(10,1:9)]), y=srg, fill=is_nrn))+
  geom_boxplot(outlier.size=.1)+
  theme_bw()+
  scale_fill_manual(values=c("grey50","red"))+
  labs(y="logcounts expr.", x="", title="HOMER1 (SRT)")

#additional plots for pde10a
vln.df = cbind.data.frame("domain"=colData(spe)[,c("domain")], "pde10a"=logcounts(spe)["PDE10A",])
load(file=here::here('plots','spatial_palette_final.rda'))

ggplot(vln.df, aes(x=domain, y=pde10a, fill=domain))+
  geom_violin(scale="width")+scale_fill_manual(values=spatial.palette)+
  theme_bw()+labs(y="Expression (logcounts)", x="", title="PDE10A")+
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

plotReducedDim(sce, dimred="UMAP", color_by="PDE10A", point_size=.1)+
  scale_color_gradient("PDE10A", low="grey90", high="black")

ggplot(mutate(as.data.frame(colData(sce)), is_nrn=broad.cell.class=="Neuron",
              srg=logcounts(sce)["PDE10A",]),
       aes(x=brnum, y=srg, fill=is_nrn))+
  geom_boxplot(outlier.size=.1)+
  theme_bw()+
  scale_fill_manual(values=c("grey50","red"))+
  labs(y="logcounts expr.", x="", title="PDE10A (snRNA-seq)")
