library(SpatialExperiment)
library(SingleCellExperiment)
library(scran)
library(dplyr)
library(ggplot2)
library(scater)
set.seed(123)

# load snRNAseq data and see whats what
load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))

load(here::here('processed-data','06_clustering','PRECAST','spe_precast_HE_domain.rda'))
levels(spe$domain)
spe$amy<-ifelse(spe$sample_id %in%
                  c('Br6423_V10B01-085_C1','Br6432_V10B01-086_B1') &
                  spe$domain %in% c('RHP'),T,F)
spe$amy<-ifelse(spe$sample_id =='Br6423_V10B01-085_A1',T,spe$amy)
spe$thal<-ifelse(spe$sample_id=='Br8325_V11A20-297_B1' & spe$domain %in% c('SUB','SUB.RHP'),T,F)

spe2 = spe[,spe$amy==FALSE & spe$thal==FALSE]

srt.degs1 <- c("KIT","COL5A2",#https://pmc.ncbi.nlm.nih.gov/articles/PMC311355/
              "POU3F1",#"MET",
              "COL24A1",#SUB
              "PART1",#SUB,
              #"IPCEF1",#RHP
              "KCNH5",#RHP
              #"ADCYAP1",#supposed to be RHP but is amygdala
              "TESPA1",#RHP
              "APOC1","SFRP2",
              "BCAS1",#WM
              "DNAH11"#Vasc/CP
              #"SLC5A5"
)

vln.df1 = cbind.data.frame("domain"=colData(spe2)[,c("domain")], as.matrix(t(logcounts(spe2)[rowData(spe2)$gene_name %in% srt.degs1,])))
colnames(vln.df1) <- c("domain",rowData(spe2)[colnames(vln.df1)[2:12],"gene_name"])

load(file=here::here('plots','spatial_palette_final.rda'))

plist1 <- lapply(srt.degs1, function(x) {
  tmp = vln.df1[,c("domain",x)]
  colnames(tmp) <- c("domain","marker")
  ggplot(tmp, aes(x=domain, y=marker, fill=domain))+
    geom_violin(scale="width")+ggtitle(x)+
    scale_fill_manual(values=spatial.palette)+
    theme_bw()+theme(axis.text.x=element_blank(),#element_text(angle=90, hjust=1), 
                     axis.title=element_blank(), 
                     legend.position="none")
}
)
do.call(gridExtra::grid.arrange, c(plist1, ncol=1))

plotGroupedHeatmap(spe2, features=srt.degs1, swap_rownames="gene_name", group="domain",
                   scale=T, center=T, 
                   cluster_rows=F, cluster_cols=F)

#srt DEGs
srt.degs <- c("KIT",#"COL5A2",#https://pmc.ncbi.nlm.nih.gov/articles/PMC311355/
              "POU3F1",#"MET",
              "COL24A1",#SUB
              "PART1",#SUB,
              #"IPCEF1",#RHP
              "KCNH5",#RHP
              #"ADCYAP1",#supposed to be RHP but is amygdala
              "TESPA1",#RHP
              #"APOC1",
              #"BCAS1",#WM
              "DNAH11"#Vasc/CP
              #"SLC5A5"
)

vln.df = cbind.data.frame("fine.cell.class"=colData(sce)[,c("fine.cell.class")], as.matrix(t(logcounts(sce)[srt.degs,])))

load(file=here::here('plots','snRNAseq_palettes.rda'))
names(sn.fine.palette) = levels(vln.df$fine.cell.class)

plist <- lapply(srt.degs, function(x) {
  tmp = vln.df[,c("fine.cell.class",x)]
  colnames(tmp) <- c("fine.cell.class","marker")
  ggplot(tmp, aes(x=fine.cell.class, y=marker, fill=fine.cell.class))+
    geom_violin(scale="width")+ggtitle(x)+
    scale_fill_manual(values=sn.fine.palette)+
    theme_bw()+theme(axis.text.x=element_blank(),#element_text(angle=90, hjust=1), 
                     axis.title=element_blank(), 
                     legend.position="none")
}
)

do.call(gridExtra::grid.arrange, c(plist, ncol=1))
