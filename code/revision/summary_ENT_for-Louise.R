library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(scater)
library(gridExtra)

set.seed(123)


####################################################################################
############################ SPATIAL REGISTRATION OF snRNAseq-DERIVED NMF PATTERNS
# code from L152-168: https://github.com/LieberInstitute/spatial_hpc/blob/main/code/revision/plots_Figure6-and-supp.R
####################################################################################

load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))
sub1 = mirrorObject(spe[,spe$sample_id=="V11U08-081_D1"], axis="v")


# visualization of NMF patterns corresponding to pyramidal snRNAseq clusters that reveal consistent laminar organization and allow for the re-annotation of clusters and more granular spatial domains 
# aka plots from Fig. 6 C and E (technically slightly different because I didn't binarize subiculum in this example)
load("plots/spatial_palette_final.rda")

broad.palette= c("Neuron"="#add294", "Neuropil"="#dfa56e", "WM"="#ff80fe", "Vasc_CSF"="#00006a")

plotVisium(sub1, spots=T, annotate="nmf27", highlight="broad.domain",
               image = FALSE, point_size=1)+coord_cartesian()+
  scale_color_manual(values=broad.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf27)),
                      labels=function(x) format(x, scientific=T, digits=1))


####################################################################################
############################ NEW ENT/RHP ANNOTATIONS BASED ON SPOT-LEVEL NMF WEIGHTS
# see extended data figures 49 C-D for context with sub, ED Fig. 50 for table of ENT/RHP, and ED Fig 51 for spot plots
# code from: https://github.com/LieberInstitute/spatial_hpc/blob/main/code/revision/plots_all-brain-spotplot_binary-nmf17.R
# also within for tables and stuff: https://github.com/LieberInstitute/spatial_hpc/blob/main/code/revision/plots_Figure6-and-supp.R
####################################################################################

# create NMF-based annotatations of ENT
#label subiculum from NMF
spe$sub1 = spe$nmf40>((max(spe$nmf40)*.95)/5)
spe$sub2 = spe$nmf54>((max(spe$nmf54)*.95)/5)
spe$subiculum = spe$sub1|spe$sub2

#label CA from NMF
spe$ca1 = spe$nmf15>((max(spe$nmf15)*.95)/5)
spe$ca3.2 = spe$nmf63>((max(spe$nmf63)*.95)/5)
spe$ca3.1 = spe$nmf11>((max(spe$nmf11)*.95)/5)
spe$ca3 = spe$ca3.1|spe$ca3.2

#label ENT from NMF
spe$ENT_sup3 = spe$nmf27>((max(spe$nmf27)*.95)/5)
spe$ENT_sup2 = spe$nmf45>((max(spe$nmf45)*.95)/5)
spe$ENT_sup1 = spe$nmf84>((max(spe$nmf84)*.95)/5)
spe$ENT_sup = spe$ENT_sup1|spe$ENT_sup2|spe$ENT_sup3
spe$ENT_L5 =  spe$nmf51>((max(spe$nmf51)*.95)/5)

#combine labels
spe$binary_label = "other"
colData(spe)[spe$subiculum==T & spe$ca1==F & spe$ca3==F & spe$ENT_sup==F & spe$ENT_L5==F,"binary_label"] = "subiculum"
colData(spe)[spe$subiculum==F & spe$ca1==T & spe$ca3==F & spe$ENT_sup==F & spe$ENT_L5==F,"binary_label"] = "CA1"
colData(spe)[spe$subiculum==F & spe$ca1==F & spe$ca3==1 & spe$ENT_sup==F & spe$ENT_L5==F,"binary_label"] = "CA3"
colData(spe)[spe$subiculum==F & spe$ca1==F & spe$ca3==F & spe$ENT_sup==T & spe$ENT_L5==F,"binary_label"] = "ENT_sup"
colData(spe)[spe$subiculum==F & spe$ca1==F & spe$ca3==F & spe$ENT_sup==F & spe$ENT_L5==T,"binary_label"] = "ENT_L5"


### use this if you want to do DE on SRT or look at expression of snRNAseq markers in corresponding spatial domains annotated from the ENT/RHP NMF patterns

####################################################################################
############################ DE RESULTS FROM snRNAseq
####################################################################################

############################ 1-vs-all: superfine cell class (60 groups)
# sn superfine DE code: https://github.com/LieberInstitute/spatial_hpc/blob/main/snRNAseq_hpc/code/revision/methods_pseudobulk-DEGs.R

reg = read.csv("snRNAseq_hpc/processed-data/revision/sn_enrichment_stats_superfine.csv", row.names=1)

tstats.df = reg[,c(grep("t_stat",colnames(reg)),237:238)] %>% 
  tidyr::pivot_longer(cols=grep("t_stat",colnames(reg), value=T), names_to="superfine.cell.class", values_to="t_stat", names_prefix="t_stat_")
fdr.df = reg[,c(grep("fdr",colnames(reg)),237:238)] %>%
  tidyr::pivot_longer(cols=grep("fdr",colnames(reg), value=T), names_to="superfine.cell.class", values_to="fdr", names_prefix="fdr_")
logfc.df = reg[,c(grep("logFC",colnames(reg)),237:238)] %>%
  tidyr::pivot_longer(cols=grep("logFC",colnames(reg), value=T), names_to="superfine.cell.class", values_to="logfc", names_prefix="logFC_")


#compile stats and limit to EXCN clusters of interest
enrich.df = left_join(tstats.df, fdr.df, by=c("superfine.cell.class","ensembl","gene")) %>%
  left_join(logfc.df, by=c("superfine.cell.class","ensembl","gene")) %>%
  filter(superfine.cell.class %in% c(#"Sub.1","Sub.2","L6.1","L2/3.1",
                                     "L6b","L6.2","L5.2","L5.1",
                                     "L2.3.2","L2.3.3","L2.3.6","L2.3.4")) %>%
  mutate(new.cell.class=factor(superfine.cell.class, levels=c(#Sub.1","Sub.2","L6.1","L2/3.1",
                                                              "L6b","L6.2","L5.2","L5.1",
                                                              "L2.3.2","L2.3.3","L2.3.6","L2.3.4"),
                               labels=c(#"Sub.1","Sub.2","Sub.3","PreS",
                                        "RHP.L6b","RHP.L6","RHP.CBLN2+","ENT.L5",
                                        "ENT.sup3","ENT.sup2b","ENT.sup2a","ENT.sup1"))) %>%
  filter(fdr<.05)

#isolating just the top 20 DE (by t stat) that are sig for only 1 group
filter(enrich.df, fdr<.05) %>% group_by(new.cell.class) %>% slice_max(n=20, t_stat) %>%
  group_by(gene) %>% add_tally(name="sig_for_n_clusters") %>%
  filter(sig_for_n_clusters==1) %>%
  group_by(new.cell.class) %>% tally(name="n_unique_sig")
# new.cell.class n_unique_sig
#RHP.L6b                  15
#RHP.L6                   15
#RHP.CBLN2+               10
#ENT.L5                   11
#ENT.sup3                 16
#ENT.sup2b                12
#ENT.sup2a                11
#ENT.sup1                 16



# ENT.sup1
filter(enrich.df, fdr<.05) %>% group_by(new.cell.class) %>% slice_max(n=20, t_stat) %>%
  group_by(gene) %>% add_tally(name="sig_for_n_clusters") %>%
  filter(sig_for_n_clusters==1, new.cell.class=="ENT.sup1") %>% 
  select(ensembl, gene, t_stat, fdr, logfc, new.cell.class)
# highlights: LINC02752, AC104689.2
# not lncRNA highlights: DSC3, CBLN1, BARX2, CAPN12, HMCN2

# ENT.sup2 NMF/spatial org comprises 2 clusters so instead of looking only for genes DE in 1 cluster, look for top genes DE in both clusters
filter(enrich.df, fdr<.05) %>% group_by(new.cell.class) %>% slice_max(n=20, t_stat) %>%
  group_by(gene) %>% add_tally(name="sig_for_n_clusters") %>%
  filter(sig_for_n_clusters==2, new.cell.class %in% c("ENT.sup2a","ENT.sup2b")) %>% 
  group_by(gene) %>% add_tally() %>% filter(n==2) %>%
  select(ensembl, gene, t_stat, fdr, logfc, new.cell.class) %>% arrange(gene)
# highlights: GLP2R, AJ009632.2

# ENT.sup3
filter(enrich.df, fdr<.05) %>% group_by(new.cell.class) %>% slice_max(n=20, t_stat) %>%
  group_by(gene) %>% add_tally(name="sig_for_n_clusters") %>%
  filter(sig_for_n_clusters==1, new.cell.class=="ENT.sup3") %>% 
  select(ensembl, gene, t_stat, fdr, logfc, new.cell.class)
# highlights: OXGR1, AVP, OR8D4, INMT, CXCL1, CTXN3, LECT2, AC092162.2, AC009081.2

#ENT L5: only lncRNA unique results in top 20 so increase to top 50
filter(enrich.df, fdr<.05) %>% group_by(new.cell.class) %>% slice_max(n=20, t_stat) %>%
  group_by(gene) %>% add_tally(name="sig_for_n_clusters") %>%
  filter(sig_for_n_clusters==1, new.cell.class=="ENT.L5") %>% 
  select(ensembl, gene, t_stat, fdr, logfc, new.cell.class) 
# top highlights: AC073578.2, AC016687.2, AC021134.1
# not lncRNA highlights in top 50: NPSR1, TYR, RORB, HTR3B


############################ 1-vs-all: subiculum/ ENT/ RHP focused
# DE run/ results saved: https://github.com/LieberInstitute/spatial_hpc/blob/main/snRNAseq_hpc/code/revision/methods_subiculum-focused-DE.R
# code for ED Fig. 52 plots: L160 start https://github.com/LieberInstitute/spatial_hpc/blob/main/snRNAseq_hpc/code/revision/plots_Figure6-and-supp.R

ent.sup = readRDS("snRNAseq_hpc/processed-data/revision/subiculum-DE_pres-superficial.rda")
sup.list =list(
  "ENT.sup1"=rownames(filter(as.data.frame(ent.sup[['ENT.sup1']]), -log10(FDR)>30 & summary.logFC>2)),
  "ENT.sup2a"=rownames(filter(as.data.frame(ent.sup[['ENT.sup2a']]), -log10(FDR)>30 & summary.logFC>2)),
  "ENT.sup2b"=rownames(filter(as.data.frame(ent.sup[['ENT.sup2b']]), -log10(FDR)>30 & summary.logFC>2)),
  "ENT.sup3"=rownames(filter(as.data.frame(ent.sup[['ENT.sup3']]), -log10(FDR)>30 & summary.logFC>2))
)
sapply(sup.list, length)
#ENT.sup1 ENT.sup2a ENT.sup2b  ENT.sup3 
#643       435       716       625
UpSetR::upset(UpSetR::fromList(sup.list))
#highlights:
### ENT.sup2a and ENT.sup2b: share 95 genes not sig in ENT.sup1 or ENT.sup3
### ENT.sup1: 220 results not sig in other superficial ENT layers
### ENT.sup3: 212 results not sig in other superficial ENT layers


############################ location of NMF feature loadings to cross reference with DE gene lists
load(file=here::here('snRNAseq_hpc','processed-data','NMF','nmf_final.rda'))
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]


####################################################################################
############################ snRNAseq visualizations 
# code from: https://github.com/LieberInstitute/spatial_hpc/blob/main/snRNAseq_hpc/code/revision/plots_Figure6-and-supp.R
####################################################################################

load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))

sce_pyr = sce[,sce$fine.cell.class %in% c("CA1/ProS","L2/3.PrS.PaS","L2/3.PrS.Ent","L5/6","L6/6b","Sub.1","Sub.2")]
sce_pyr$new.cell.class = factor(sce_pyr$superfine.cell.class,
                                levels=c("L2/3.5","CA1","ProS",
                                         "Sub.1","Sub.2","L6.1","L2/3.1",
                                         "L6b","L6.2","L5.2","L5.1",
                                         "L2/3.2","L2/3.3","L2/3.6","L2/3.4"),
                                labels=c("L2/3.5","CA1","ProS",
                                         "Sub.1","Sub.2","Sub.3","PreS",
                                         "RHP.L6b","RHP.L6","RHP.CBLN2+","ENT.L5",
                                         "ENT.sup3","ENT.sup2b","ENT.sup2a","ENT.sup1")
)


#dot plots
plotDots(sce_pyr, features=c("RORB", "AC021134.1",
                              "AC092162.2", "CTXN3", "AC009081.2",
                              "GLP2R", "AJ009632.2",
                             #"LINC02752", # great for ENT.sup1 but too high of expr for rest of these genes (changes limits)
                             "AC104689.2", "BARX2", "HMCN2"), group="new.cell.class")+
  scale_color_viridis_c(option="F", direction=-1)+
  labs(size="prop. nuclei", color="Avg. expr.")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), text=element_text(size=11))


#TSNE (identical to that in Fig. 6 D and F)
set.seed(123)
reducedDim(sce_pyr, "TSNE") = calculateTSNE(sce_pyr, dimred="MNN", n_dimred=50)

#nmfs
nmf.list = c("nmf15","nmf32","nmf40","nmf54","nmf84","nmf45","nmf27","nmf51","nmf68","nmf22","nmf53","nmf65")

plist1 = lapply(nmf.list, function(x) 
  plotReducedDim(sce_pyr, dimred="TSNE", color_by=x, point_size=.1)+
    scale_color_viridis_c(option="F", direction=-1)+ggtitle(x)+
    theme(legend.position="none", axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(),
          plot.title=element_text(hjust=.5, size=14))
)
do.call(gridExtra::grid.arrange, c(plist1, ncol=4))

#marker genes
plist2 = lapply(c("RORB", "AC021134.1",
                 "AC092162.2", "CTXN3", "AC009081.2",
                 "GLP2R", "AJ009632.2",
                 #"LINC02752", # great for ENT.sup1 but too high of expr for rest of these genes (changes limits)
                 "AC104689.2", "BARX2", "HMCN2"), function(x) 
  plotReducedDim(sce_pyr, dimred="TSNE", color_by=x, point_size=.1)+
    scale_color_viridis_c(option="F", direction=-1)+ggtitle(x)+
    theme(legend.position="none", axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(),
          plot.title=element_text(hjust=.5, size=14))
)
do.call(gridExtra::grid.arrange, c(plist2, ncol=3))
