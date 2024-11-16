library(SpatialExperiment)
library(scater)
library(dplyr)
library(ggplot2)
set.seed(123)

#umap for schematic
load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))
unique(sce$brnum)
plotReducedDim(sce, dimred="UMAP", color_by="sort", point_size=.3)+
  scale_color_manual(values=c("#875faa","#64a65a"))

plotReducedDim(sce[,sce$sort=="PI+NeuN+"], dimred="UMAP", color_by="sort", point_size=.5)+
  scale_color_manual(values=c("#64a65a"))

################ spot plots and manual annotation markers
load(here::here('processed-data','06_clustering','PRECAST','spe_precast_HE_domain.rda'))
table(spe$VSPG)

#spot plot for MBP
sub = spe[,spe$sample_id %in% c("Br3942_V11L05-333_C1","Br3942_V11L05-333_D1")]
plotSpots(sub, annotate=rownames(spe)[rowData(spe)$gene_name=="MBP"],
          sample_id="sample_id", assay="logcounts")+
  scale_color_gradient(low="grey90", high="black")+
  theme_void()+theme(strip.text=element_blank())


spe$amy<-ifelse(spe$sample_id %in%
                  c('Br6423_V10B01-085_C1','Br6432_V10B01-086_B1') &
                  spe$domain %in% c('RHP'),T,F)
spe$amy<-ifelse(spe$sample_id =='Br6423_V10B01-085_A1',T,spe$amy)
spe$thal<-ifelse(spe$sample_id=='Br8325_V11A20-297_B1' & spe$domain %in% c('SUB','SUB.RHP'),T,F)

spe$cluster = as.character(spe$cluster)
spe$cluster = ifelse(spe$amy==TRUE, "Amygdala",spe$cluster)
spe$cluster = ifelse(spe$thal==TRUE, "Thalamus",spe$cluster)

#i don't need to include the marker genes cause those are in Fig 2D
genes = c(#"SHOX2",
          "TCF7L2",
          "GAD2","SST",
          #"SLC17A6",
          #"CDH22","OPRM1","CACNG4",
          "RORB","CUX2",#"SATB2","TLE4",
          "NTS","FN1",
          "MPPED1","FIBCD1",
          "TSPAN18","NECTIN3","AMPH",
          #"SFRP2",
          "PROX1","SEMA5A",
          #"SOX2","SLC1A2",
          "MT-ND5","MT-ATP8","GFAP",
          "MOG","MOBP","PLP1","PRLR","TTR"#,"PDGFRB","TPM2"
          )

spe$cluster_f = factor(spe$cluster, levels=c("Thalamus","GABA","Amygdala",
                                             "RHP","SUB.RHP","SUB",
                                             "CA1.1","CA1.2","CA2.4.1","CA2.4.2",
                                             "GCL",
                                             "ML","SL.SR","SR.SLM","SLM.SGZ",
                                             "WM.1","WM.2","WM.3","Vascular","Choroid"))

plotGroupedHeatmap(spe, features=genes, group="cluster_f", swap_rownames="gene_name",
                   center=T, scale=T, zlim=c(-2,4),
                   cluster_rows=F, cluster_cols=F,
                   #colour=viridisLite::magma(n=10))
                   colour=viridisLite::magma(n=30)[c(rep(1,5),2:30)],
                   angle_col=90)

load(file = here::here("processed-data","nnSVG","nnSVG_gene_lists.rda"))

intersect(genes, rowData(spe)[nnSVG,"gene_name"])
setdiff(genes, rowData(spe)[nnSVG,"gene_name"])


tmp = spe[,spe$sample_id=="Br6471_V11L05-335_C1"]
tmp = spe[,spe$sample_id=="Br3942_V11L05-333_D1"]
library(ggspavis)
p1 <- plotSpots(tmp, annotate=rownames(spe)[rowData(spe)$gene_name=="MOBP"], assay_name="logcounts",
                pal=c("grey90","black"))+
  theme_void()+ggtitle("MOBP")
p2 <- plotSpots(tmp, annotate=rownames(spe)[rowData(spe)$gene_name=="GFAP"], assay_name="logcounts",
                pal=c("grey90","black"))+
  theme_void()+ggtitle("GFAP")
p3 <- plotSpots(tmp, annotate=rownames(spe)[rowData(spe)$gene_name=="AMPH"], assay_name="logcounts",
                pal=c("grey90","black"))+
  theme_void()+ggtitle("AMPH")
p4 <- plotSpots(tmp, annotate=rownames(spe)[rowData(spe)$gene_name=="SEMA5A"], assay_name="logcounts",
                pal=c("grey90","black"))+
  theme_void()+ggtitle("SEMA5A")+scale_color_gradient(low="grey90", high="black", breaks=c(0,1,2))

gridExtra::grid.arrange(p2, p1, p3, p4, ncol=2)
