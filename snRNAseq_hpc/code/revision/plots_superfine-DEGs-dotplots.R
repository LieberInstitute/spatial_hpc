library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scater)
set.seed(123)

load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))

reg = read.csv("snRNAseq_hpc/processed-data/revision/sn_enrichment_stats_superfine.csv", row.names=1)

tstats.df = reg[,c(grep("t_stat",colnames(reg)),237:238)] %>% 
  tidyr::pivot_longer(cols=grep("t_stat",colnames(reg), value=T), names_to="superfine.cell.class", values_to="t_stat", names_prefix="t_stat_")
fdr.df = reg[,c(grep("fdr",colnames(reg)),237:238)] %>%
  tidyr::pivot_longer(cols=grep("fdr",colnames(reg), value=T), names_to="superfine.cell.class", values_to="fdr", names_prefix="fdr_")
logfc.df = reg[,c(grep("logFC",colnames(reg)),237:238)] %>%
  tidyr::pivot_longer(cols=grep("logFC",colnames(reg), value=T), names_to="superfine.cell.class", values_to="logfc", names_prefix="logFC_")

enrich.df = left_join(tstats.df, fdr.df, by=c("superfine.cell.class","ensembl","gene")) %>%
  left_join(logfc.df, by=c("superfine.cell.class","ensembl","gene"))
#Macro/Tcell" "PC/SMC
enrich.df = left_join(enrich.df, 
                      distinct(as.data.frame(colData(sce)[,c("fine.cell.class","superfine.cell.class")])) %>%
                        mutate(superfine.cell.class=gsub("/",".",superfine.cell.class)),
                      by="superfine.cell.class")
namecode = distinct(enrich.df, fine.cell.class, superfine.cell.class) %>% mutate(new_name=paste(fine.cell.class, superfine.cell.class, sep="_"))

up.df = filter(enrich.df, fdr<.0001 & logfc>2)
distinct(up.df, gene, superfine.cell.class) %>% group_by(gene) %>% add_tally() %>%
  filter(n==1) %>% group_by(superfine.cell.class) %>% tally()

hpc = sce[,sce$fine.cell.class %in% c("CA1/ProS","CA2-4","Sub.1","Sub.2","MC")]
hpc$superfine.cell.class = droplevels(hpc$superfine.cell.class)
hpc.genes <- c("FBLN5","CARTPT", #"MORC1", #MC
               "TSPAN18", #CA3
               "TRPS1","NUP93",
               "COL5A2",#"MCTP2",
               "CHST8",
               "FN1",#Sub1
               "CNTN6", "RPRM")
hpc.plot <- plotDots(hpc, features=hpc.genes, group="superfine.cell.class", assay.type = "logcounts")+
  scale_color_gradient(low="grey90", high="black")+ggtitle("HPC")+
  scale_size(range=c(0,3))+labs(size="prop.",color="avg. expr.")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), axis.title=element_blank(),
        legend.text = element_text(size = 6), legend.key.size = unit(10,"pt"),
        legend.title = element_text(size=8))

gaba.big = sce[,sce$fine.cell.class %in% c("GABA.MGE","GABA.CGE","GABA.LAMP5","GABA.PENK","Cajal")]
gaba.big$superfine.cell.class = factor(gaba.big$superfine.cell.class, 
                                   levels=c("PENK","CORT","SST","PV.FS","C1QL1","CRABP1",
                                            "LAMP5.MGE","LAMP5.CGE","CXCL14",
                                            "HTR3A","VIP","Cajal"
                                   ))
gaba.big.genes = c("GAD1","GAD2","NPY","LHX6","NXPH1","LAMP5","ADARB2","CCK","RELN","TP73")
gaba.big.plot <- plotDots(gaba.big, features=gaba.big.genes, group="superfine.cell.class")+
  scale_color_gradient(low="grey90", high="black")+ggtitle("GABA (high expr)")+
  scale_size(range=c(0,3))+labs(size="prop.",color="avg. expr.")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), axis.title=element_blank(),
        legend.text = element_text(size = 6), legend.key.size = unit(10,"pt"),
        legend.title = element_text(size=8))



gaba.small = sce[,sce$fine.cell.class %in% c("GABA.MGE","GABA.CGE","GABA.LAMP5","GABA.PENK")]
gaba.small$superfine.cell.class = factor(gaba.small$superfine.cell.class, 
                                   levels=c("PENK","CORT","SST","PV.FS","C1QL1","CRABP1",
                                            "LAMP5.MGE","LAMP5.CGE","CXCL14",
                                            "HTR3A","VIP"))
gaba.small.genes = c("PENK","MEIS2","CORT","SST",
                     "TAC1","PVALB","C1QL1",
                     "SFTA3","CRABP1",#"ID2",
                     "CCK",
                     "HAPLN1","SV2C",#"CPLX3",
                     "KIT","CNR1",
                     "CALB2","VIP")
gaba.small.plot <- plotDots(gaba.small, features=gaba.small.genes, group="superfine.cell.class")+
  scale_color_gradient(low="grey90", high="black")+ggtitle("GABA (lower expr)")+
  scale_size(range=c(0,3))+labs(size="prop.",color="avg. expr.")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), axis.title=element_blank(),
        legend.text = element_text(size = 6), legend.key.size = unit(10,"pt"),
        legend.title = element_text(size=8))



cort = sce[,sce$fine.cell.class %in% c("L2/3.PrS.PaS","L2/3.PrS.Ent","L5/6","L6/6b")]
cort$superfine.cell.class = factor(cort$superfine.cell.class, 
                                   levels=c("L2/3.2","L2/3.4","L2/3.3","L2/3.6",
                                            "L2/3.1","L2/3.5",
                                            "L5.1","L5.2","L6.1","L6b","L6.2"))
cort.genes = c("CBLN4","TRPS1","CUX2","NTNG1","TSHZ2","ABCA12","VIPR2","ADCYAP1","SH3RF2",
               "SATB2","RORB","THEMIS","CBLN2","NXPH3","EFHD2","SEMA3E","TLE4","FOXP2")
cort.plot <- plotDots(cort, features=cort.genes, group="superfine.cell.class")+
  scale_color_gradient(low="grey90", high="black")+ggtitle("RHP")+
  scale_size(range=c(0,3))+labs(size="prop.",color="avg. expr.")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), axis.title=element_blank(),
        legend.text = element_text(size = 6), legend.key.size = unit(10,"pt"),
        legend.title = element_text(size=8))




amy = sce[,sce$fine.cell.class %in% c("Amy","L2/3.PrS.Ent","L2/3.PrS.PaS","HATA","L6/6b","L5/6")]
amy$superfine.cell.class <- droplevels(amy$superfine.cell.class)
amy$new_name = factor(amy$superfine.cell.class, levels=levels(amy$superfine.cell.class),
                      labels=c("L6/6b","L6/6b","L6/6b","L5/6","L5/6",
                               "L2/3.2","L2/3.2","L2/3.1",
                               "L2/3.2","L2/3.2","L2/3.1",
                               "HATA","AHi.1","AHi.2","AHi.3","AHi.4"))
amy$new_name = factor(as.character(amy$new_name),
                      levels=c("HATA","AHi.1","AHi.2","AHi.3","AHi.4",
                               "L6/6b","L5/6","L2/3.2","L2/3.1"))
amy.genes = c("POSTN","ESR1","PTGER3",
  "TACR3",
  "FREM3",
  "DCSTAMP",
  "NPFFR2","PLEKHG4B",
  "CORIN","PAPPA2","PSME4")
amy.plot <- plotDots(amy, features=amy.genes, group="new_name")+
  scale_color_gradient(low="grey90", high="black")+ggtitle("Amygdala")+
  scale_size(range=c(0,3))+labs(size="prop.",color="avg. expr.")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), axis.title=element_blank(),
        legend.text = element_text(size = 6), legend.key.size = unit(10,"pt"),
        legend.title = element_text(size=8))



gcs = sce[,sce$fine.cell.class=="GC"]
gcs$superfine.cell.class <- factor(gcs$superfine.cell.class, levels=c("GC.3","GC.2","GC.1","GC.5","GC.4"))
gc.genes = c("ACVR2A","CNGB1","COL25A1",
  "ACVR1","ITGA8",
  "ARHGEF3","BDNF",
  "ACVR1C"
)

gcs.plot <- plotDots(gcs, features=gc.genes, group="superfine.cell.class")+
  scale_color_gradient(low="grey90", high="black")+ggtitle("Granule cells")+
  scale_size(range=c(0,3))+labs(size="prop.",color="avg. expr.")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), axis.title=element_blank(),
        legend.text = element_text(size = 6), legend.key.size = unit(10,"pt"),
        legend.title = element_text(size=8))



glia <- sce[,sce$superfine.cell.class %in% c("Astro.1","Astro.2","Astro.3","Oligo.1","Oligo.2","OPC","COP",
                                             "Micro.1","Micro.2","Macro/Tcell","Ependy",
                                             "CP.1","CP.2","CP.3","Endo.2","Endo.1","PC/SMC","VLMC")]
glia$superfine.cell.class = factor(glia$superfine.cell.class,
                       levels=c("Astro.2","Astro.3","Astro.1","Oligo.1","Oligo.2","OPC","COP",
                                "Micro.1","Micro.2","Macro/Tcell","Ependy",
                                "CP.1","CP.2","CP.3","Endo.2","Endo.1","PC/SMC","VLMC"))
glia.genes= c("TNC","DOCK7","ETNPPL","RORB",
              "PLP1","MOBP","MOG","BCAS1","LHFPL3",
              "GPR17",#COPs = differentiation committed, GPR17 from this paper: https://www.hznu.edu.cn/upload/resources/file/2023/09/18/7791590.pdf
              "P2RY12","APBB1IP",
              "CFAP299","DNAH11",
              #"OTX2-AS1",
              "SLC13A4",
              "IGFBP7",
              "DKK2","MECOM","UTRN","ABCB1",#endo #ABCB1: This protein also functions as a transporter in the blood-brain barrier
              "DLC1","LAMA2","CEMIP")#pc/smc, vlmc

glia.plot <- plotDots(glia, features=glia.genes, group="superfine.cell.class")+
  scale_color_gradient(low="grey90", high="black")+ggtitle("Glia")+
  scale_size(range=c(0,3))+labs(size="prop.",color="avg. expr.")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), axis.title=element_blank(),
        legend.text = element_text(size = 6), legend.key.size = unit(10,"pt"),
        legend.title = element_text(size=8))





#glia
c(length(glia.genes),length(levels(glia$superfine.cell.class))) #23,18

#amy
c(length(amy.genes),length(levels(amy$new_name))) #11,9
#hpc
c(length(hpc.genes),length(levels(hpc$superfine.cell.class))) #10,8


#glia 23x18
#amy (11,9) + hpc (10,8)
mat1 = matrix(1, nrow=23, ncol=18)
mat2 = cbind(rbind(matrix(2, nrow=10, ncol=8), rep(NA, 8)), cbind(rep(NA,11), matrix(3,nrow=11, ncol=9)))

plot1 = glia.plot
plot2 = hpc.plot
plot3 = amy.plot

left.mat = rbind(mat1, rep(NA, 18), mat2)


#cort
c(length(cort.genes),length(levels(cort$superfine.cell.class))) #18,11
#gaba small
c(length(gaba.small.genes),length(levels(gaba.small$superfine.cell.class))) #16,11

#next to that is cort on top of gaba small (34 tall, both 11 wide)
mat3 = rbind(matrix(4, nrow=18, ncol=11), rep(NA, 11), matrix(5,nrow=16, ncol=11))

plot4 = cort.plot
plot5 = gaba.small.plot

top.mat = cbind(left.mat, rep(NA, 35), mat3)

#gaba big
c(length(gaba.big.genes),length(levels(gaba.big$superfine.cell.class))) #10,12
#gcs
c(length(gc.genes),length(levels(gcs$superfine.cell.class))) #8,5

bot.left = cbind(matrix(7, nrow=10, ncol=12), rep(NA,10), 
                 rbind(rep(NA, 10), matrix(8, nrow=8, ncol=10), rep(NA, 10)))

plot6 = gaba.big.plot
plot7 = gcs.plot

dim(top.mat) #35, 30
#dim(bot.left) #10, 18
dim(bot.left) #extended GCs: 10,23

bot.right = matrix(NA, nrow=10, ncol=30-23)

bot.mat = cbind(bot.left, bot.right)

laymat = rbind(top.mat, rep(NA, ncol(top.mat)), bot.mat)

gridExtra::grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, layout_matrix=laymat)

dim(laymat)

pdf(file = "snRNAseq_hpc/plots/revision/dotplot_superfine_DEGs.pdf",
    width=10, height=15)
gridExtra::grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, layout_matrix=laymat)
dev.off()
