library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scater)
set.seed(123)

load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))
unique(sce$brnum)
#adjust levels for plot order
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

###cell type group dotplots

#HPC
#lots of filtering to find genes for each cluster that were significant and worked well together when plotted (same range/ scale) 
filter(up.df, superfine.cell.class=="ProS") %>% arrange(desc(t_stat)) %>% pull(gene)

hpc = sce[,sce$fine.cell.class %in% c("CA1/ProS","CA2-4","Sub.1","Sub.2","MC")]
plotGroupedHeatmap(hpc, features=c("FIBCD1","FNDC1","MCTP2","DMRTA1","MORC1","GUCA1C","GHSR",
                                   "FRMD7","MLPH","TPH2","GABRQ","RGS14","SCGN","NR3C2","RPRM","CNTN6",
                                   "ABCC12","CHRM5","ITGA9","SLC6A5","OSR1","KREMEN2"),
                   group="superfine.cell.class",
                   center=T, scale=T)

plotDots(hpc, features=c("FBLN5","CARTPT", #"MORC1", #MC
                         "TSPAN18", #CA3
                         "TRPS1",#"NECTIN3",#"HPCA",#"DOC2B","GABRQ",
                         #"ARHGEF3",
                         "NUP93","PRLR",
                         #"OLFML2B","RGS14", "SCGN", #CA2
                         # "FIBCD1", "MCTP2", #CA1
                         "COL5A2",
                          "FNDC1","POU3F1",#ProS
                         "CHST8",
                         # "CHRM5",
                        "COL24A1", # subiculum, human snrnaseq https://www.cell.com/neuron/fulltext/S0896-6273(21)00329-9 figure S7
                         "FN1",#Sub1
                          "CNTN6", "RPRM"), #Sub.2
                   group="superfine.cell.class")+
  scale_color_viridis_c(option="F", direction=-1)

#GABA
check.genes = filter(up.df, superfine.cell.class %in% c("VIP","HTR3A")) %>% group_by(gene) %>% tally() %>% filter(n==2) %>% pull(gene)
intersect(filter(up.df, gene %in% check.genes & superfine.cell.class=="HTR3A") %>% slice_max(n=20, t_stat) %>% pull(gene),
          filter(up.df, gene %in% check.genes & superfine.cell.class=="VIP") %>% slice_max(n=20, t_stat)%>% pull(gene))
#big scale
gaba = sce[,sce$fine.cell.class %in% c("GABA.MGE","GABA.CGE","GABA.LAMP5","GABA.PENK","Cajal")]
gaba$superfine.cell.class = factor(gaba$superfine.cell.class, 
                                   levels=c("PENK","CORT","SST","PV.FS","C1QL1","CRABP1",
                                            "LAMP5.MGE","LAMP5.CGE","CXCL14",
                                            "HTR3A","VIP","Cajal"
                                            ))
plotDots(gaba, features=c("GAD1","GAD2","NPY","LHX6","NXPH1",
                          "LAMP5","ADARB2","CCK","RELN","TP73"),
         group="superfine.cell.class")+
  scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=45, hjust=1))
#small scale
gaba = sce[,sce$fine.cell.class %in% c("GABA.MGE","GABA.CGE","GABA.LAMP5","GABA.PENK")]
gaba$superfine.cell.class = factor(gaba$superfine.cell.class, 
                                   levels=c("PENK","CORT","SST","PV.FS","C1QL1","CRABP1",
                                            "LAMP5.MGE","LAMP5.CGE","CXCL14",
                                            "HTR3A","VIP"))
plotDots(gaba, features=c("PENK",#"MEIS2",
                          "CORT","SST",
                          "TAC1",
                          "PVALB","C1QL1",
                          "SFTA3","CRABP1",#"ID2",
                          "CCK",
                          "HAPLN1","SV2C",#"CPLX3",
                          "KIT","CNR1",
                          "CALB2","VIP"),
         group="superfine.cell.class")+
  scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=45, hjust=1))

#cortical
filter(up.df, fine.cell.class %in% c("L2/3.PrS.PaS","L2/3.PrS.Ent")) %>% slice_max(n=50, t_stat) %>% pull(gene)

filter(up.df, fine.cell.class %in% c("L2/3.PrS.PaS")) %>% group_by(gene) %>% tally() %>% filter(n>1) %>% pull(gene)

cort = sce[,sce$fine.cell.class %in% c("L2/3.PrS.PaS","L2/3.PrS.Ent","L5/6","L6/6b")]
cort$new_name = paste(cort$fine.cell.class, cort$superfine.cell.class, sep="_")
cort$new_name = factor(cort$new_name, levels=c("L2/3.PrS.Ent_L2/3.2","L2/3.PrS.Ent_L2/3.4","L2/3.PrS.Ent_L2/3.3","L2/3.PrS.Ent_L2/3.6",
                                               "L2/3.PrS.PaS_L2/3.1","L2/3.PrS.PaS_L2/3.5",
                                               "L5/6_L5.1","L5/6_L5.2","L6/6b_L6.1","L6/6b_L6b","L6/6b_L6.2"))
plotDots(cort, features=c("ADCYAP1","PIK3C2G","GPR88","CBLN4","NRAD1","CHRNA3","GPR6","GLP2R","IQCJ",
                          "EMX1","HTR2A","KCNH4","KCNV1","KLF5","SH3RF2"),
         group="new_name")+
  scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=45, hjust=1))

filter(up.df, fine.cell.class %in% c("L5/6")) %>% group_by(gene) %>% tally() %>% filter(n>1) %>% pull(gene)


plotDots(cort, features=c("CBLN4","TRPS1","CUX2","NTNG1",
                          "TSHZ2",
                          "ABCA12","VIPR2",#"MEF2C","GRASP","ZMAT4",
                          "ADCYAP1","SH3RF2",
                          "SATB2",
                          "RORB","THEMIS","CBLN2",
                          "NXPH3","EFHD2",
                          "SEMA3E","TLE4","FOXP2","COL24A1"),
         group="new_name")+
  scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=45, hjust=1))

##HATA and AHi
amy = sce[,sce$fine.cell.class %in% c("Amy","L2/3.PrS.Ent","L2/3.PrS.PaS","HATA","L6/6b","L5/6")]
amy$superfine.cell.class <- droplevels(amy$superfine.cell.class)
amy$new_name = factor(amy$superfine.cell.class, levels=levels(amy$superfine.cell.class),
                      labels=c("L6/6b","L6/6b","L6/6b","L5/6","L5/6",
                               "L2/3.PrS.Ent","L2/3.PrS.Ent","L2/3.PrS.PaS",
                               "L2/3.PrS.Ent","L2/3.PrS.Ent","L2/3.PrS.PaS",
                               "HATA","AHi.1","AHi.2","AHi.3","AHi.4"))
amy$new_name = factor(as.character(amy$new_name),
                      levels=c("HATA","AHi.1","AHi.2","AHi.3","AHi.4",
                               "L6/6b","L5/6","L2/3.PrS.Ent","L2/3.PrS.PaS"))
plotDots(amy, features=c(#"SLC17A6","SCN11A","M1AP",
                         "POSTN","ESR1","PTGER3",
                         "TACR3",
                         "FREM3",
                         "DCSTAMP",
                         "NPFFR2","PLEKHG4B",
                         "CORIN","PAPPA2","PSME4"),
         group="new_name")+
  #scale_color_viridis_c(option="F", direction=-1)+
  scale_color_gradientn(colors=viridisLite::rocket(6)[6:1], values=c(seq(0,.7,length.out=5),1))+
  theme(axis.text.x=element_text(angle=45, hjust=1))

#GCs 
gcs = sce[,sce$fine.cell.class=="GC"]
gcs$superfine.cell.class <- droplevels(gcs$superfine.cell.class)
gcs$new_name = factor(gcs$superfine.cell.class, levels=c("GC.3","GC.2","GC.1","GC.5","GC.4"))

  filter(up.df, superfine.cell.class %in% c("GC.1","GC.2","GC.3","GC.5","GC.4")) %>% 
  group_by(gene) %>% add_tally() %>% filter(superfine.cell.class=="GC.2") %>% 
  arrange(desc(logfc)) %>% pull(gene)

plotDots(gcs, features=c(#"AKAIN1","PMEPA1",
                         "ACVR2A","CNGB1","COL25A1",
                         "ACVR1","ITGA8",
                         "ARHGEF3","BDNF",
                         "ACVR1C"
                         ),
  group="new_name")+
  scale_color_viridis_c(option="F", direction=-1)

## glia
glia <- sce[,sce$superfine.cell.class %in% c("Astro.1","Astro.2","Astro.3","Oligo.1","Oligo.2","OPC","COP",
                                             "Micro.1","Micro.2","Macro/Tcell","Ependy",
                                             "CP.1","CP.2","CP.3","Endo.2","Endo.1","PC/SMC","VLMC")]
glia$superfine.cell.class <- droplevels(glia$superfine.cell.class)
glia$new_name = factor(glia$superfine.cell.class,
                       levels=c("Astro.2","Astro.3","Astro.1","Oligo.1","Oligo.2","OPC","COP",
                                "Micro.1","Micro.2","Macro/Tcell","Ependy",
                                "CP.1","CP.2","CP.3","Endo.2","Endo.1","PC/SMC","VLMC"))

glia.genes = filter(up.df, superfine.cell.class %in% glia$superfine.cell.class) %>% 
  group_by(superfine.cell.class) %>% slice_max(n=5, t_stat) %>% pull(gene) %>% unique()

#pericyte/smooth muscle cell: RGS5, CNN1, CD73, CD146, CD90 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6428685/
#VLMC = vascular leptomeningeal cells: Col1A1, col1A2, LUM, DCN, FBLN1 https://www.nature.com/articles/s41586-018-0654-5 
#
setdiff(filter(up.df, superfine.cell.class=="Endo.1") %>% pull(gene),filter(up.df, superfine.cell.class=="Endo.2") %>% pull(gene))

plotDots(glia, features=c("TNC","DOCK7","ETNPPL","RORB",
                          "PLP1","MOBP","MOG","BCAS1","LHFPL3",
                          "GPR17",#COPs = differentiation committed, GPR17 from this paper: https://www.hznu.edu.cn/upload/resources/file/2023/09/18/7791590.pdf
                          "P2RY12","APBB1IP",
                          "CFAP299","DNAH11",
                          #"OTX2-AS1",
                          "SLC13A4",
                          "IGFBP7",
                          "DKK2","MECOM","UTRN","ABCB1",#endo #ABCB1: This protein also functions as a transporter in the blood-brain barrier
                          "DLC1","LAMA2","CEMIP"),#pc/smc, vlmc
         group="new_name")+
  scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=45, hjust=1))
