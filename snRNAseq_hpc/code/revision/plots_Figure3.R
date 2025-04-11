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
#              "DNAH11"#Vasc/CP
              #"SLC5A5"
)

#snRNAseq violin plots
sce$fine.cell.class = factor(sce$fine.cell.class, levels=levels(sce$fine.cell.class),
                             labels=c(levels(sce$fine.cell.class)[1:7], "L2/3.1","L2/3.2",
                                      levels(sce$fine.cell.class)[10:24]))

vln.df = cbind.data.frame("fine.cell.class"=colData(sce)[,c("fine.cell.class")], as.matrix(t(logcounts(sce)[srt.degs,])))
long.df = tidyr::pivot_longer(vln.df, all_of(srt.degs), names_to="marker.gene", values_to="expr") %>%
  mutate(marker.gene=factor(marker.gene, levels=srt.degs))

load(file=here::here('plots','snRNAseq_palettes.rda'))
names(sn.fine.palette) = levels(sce$fine.cell.class)

p1 <- ggplot(long.df, aes(x=fine.cell.class, y=expr, fill=fine.cell.class))+
  geom_violin(scale="width")+scale_fill_manual(values=sn.fine.palette)+
  facet_wrap(vars(marker.gene), ncol=1, scale="free_y")+
  theme_bw()+labs(y="Expression (logcounts)", x="")+
  theme(legend.position="none", text=element_text(size=12, color="black"), panel.grid=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=.5), axis.title.x=element_blank(),
        strip.text=element_text(size=12),
        strip.background=element_rect(fill="white", color="transparent"))

pdf(file = "snRNAseq_hpc/plots/revision/Figure3_sn-srt-degs_no-DNAH11.pdf",
    width=5.5, height=7.7)
p1
dev.off()

#all superfine cell classes heatmap
sce$superfine.cell.class = factor(sce$superfine.cell.class,
                                  levels=c("GC.3","GC.1","GC.2","GC.4","GC.5",
                                           "MC","CA3.1","CA3.2","CA2","CA1","ProS","Sub.1","Sub.2",
                                           "L6.1","L6.2","L6b","L5.2","L5.1",
                                           "L2/3.3","L2/3.6","L2/3.4","L2/3.2",
                                           "L2/3.1","L2/3.5",
                                           "HATA","AHi.1","AHi.2","AHi.3","AHi.4","Thal",
                                           "Cajal","CXCL14","HTR3A","VIP",
                                           "LAMP5.CGE","LAMP5.MGE","CRABP1",
                                           "C1QL1","PV.FS","SST","CORT","PENK",
                                           "Astro.1","Astro.2","Astro.3",
                                           "Oligo.1","Oligo.2","OPC","COP","Micro.1","Micro.2","Macro/Tcell",
                                           "Ependy","CP.1","CP.2","CP.3","Endo.2","Endo.1","PC/SMC","VLMC"))
hm = plotGroupedHeatmap(sce, features=c("ACVR2A","PROX1","ACVR1","ACVR1C","BDNF",
                                   "CARTPT","TRPS1","TSPAN18","KIT",#"NUP93",
                                   "FIBCD1","COL5A2","CHRM5","CHST8",
                                   "FN1","COL24A1","CNTN6",
                                   "NXPH3","EFHD2","TLE4","SATB2",
                                   "CBLN2","RORB",
                                   "CBLN4","NTNG1",
                                   "TESPA1","TSHZ2","VIPR2","SH3RF2",
                                   "ESR1","PTGER3","NPFFR2","PAPPA2","TCF7L2",
                                   "TP73","RELN",
                                   "ADARB2","VIP","LAMP5","C1QL1","NXPH1","LHX6","TAC1","NPY","PENK",
                                   "DOCK7","TNC",
                                   "PLP1","BCAS1","LHFPL3","GPR17",
                                   "P2RY12","APBB1IP",
                                   "CFAP299","DNAH11","SLC13A4",
                                   "IGFBP7","MECOM","DLC1","LAMA2"
                                   ),
                   group="superfine.cell.class",
                   scale=T, center=T,
                   cluster_rows=F, cluster_cols=F, zlim=c(-2,4),
                   colour=c("#FFF",RColorBrewer::brewer.pal(9,"Greys")),#viridisLite::magma(n=30)[c(rep(1,5),2:30)],
                   angle_col=90, silent=T) 
ggsave("snRNAseq_hpc/plots/revision/Figure3_heatmap.pdf", hm, bg="white", width=16, height=10, units="in")
