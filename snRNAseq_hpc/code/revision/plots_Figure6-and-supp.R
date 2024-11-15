library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scater)

set.seed(123)

load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))

## tsne
sce_pyr = sce[,sce$fine.cell.class %in% c("CA1/ProS","L2/3.PrS.PaS","L2/3.PrS.Ent","L5/6","L6/6b","Sub.1","Sub.2")]
reducedDim(sce_pyr, "TSNE") = calculateTSNE(sce_pyr, dimred="MNN", n_dimred=50)

nmf.list = c("nmf15","nmf32","nmf40","nmf54","nmf84","nmf45","nmf27","nmf51","nmf68","nmf22","nmf53","nmf65")

plist = lapply(nmf.list, function(x) 
  plotReducedDim(sce_pyr, dimred="TSNE", color_by=x, point_size=.3)+
  scale_color_viridis_c(option="F", direction=-1)+ggtitle(x)+
  theme(legend.position="none", axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(),
        plot.title=element_text(hjust=.5, size=18))
)
do.call(gridExtra::grid.arrange, c(plist, ncol=4))

##############################
######### L6.1 as part of subiculum
##############################

sce_subset = sce[,sce$mid.cell.class=="ExcN"]
sce_subset$superfine.cell.class = droplevels(sce_subset$superfine.cell.class)

ggplot(filter(as.data.frame(colData(sce_subset)), superfine.cell.class=="L6b"),
       aes(x=nmf53, y=nmf65, color=nmf65>.00035))+
  geom_point()+theme_bw()+ggtitle("L6b only")
tmp = sce_subset[,sce_subset$superfine.cell.class=="L6b"]
tmp$nmf65pos = tmp$nmf65>.00035


colData(sce_subset)$is_deep = ifelse(sce_subset$fine.cell.class %in% c("L5/6","L6/6b","Sub.1","Sub.2"), sce_subset$new.classes, "other")
sce_deep = sce_subset[,sce_subset$is_deep!="other"]
sce_deep$new.classes = as.character(sce_deep$superfine.cell.class)
colData(sce_deep)[colnames(tmp)[tmp$nmf65pos],"new.classes"] = "L6b_nmf65"
sce_deep$new.classes = factor(sce_deep$new.classes, levels=c("Sub.1","Sub.2","L6.1","L6b_nmf65","L6b","L6.2","L5.1","L5.2"))

#load nmf patters
load(file=here::here('snRNAseq_hpc','processed-data','NMF','nmf_final.rda'))
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

#load sn DEGs
reg = read.csv("snRNAseq_hpc/processed-data/revision/sn_enrichment_stats_superfine.csv", row.names=1)

tstats.df = reg[,c(grep("t_stat",colnames(reg)),237:238)] %>% 
  tidyr::pivot_longer(cols=grep("t_stat",colnames(reg), value=T), names_to="superfine.cell.class", values_to="t_stat", names_prefix="t_stat_")
fdr.df = reg[,c(grep("fdr",colnames(reg)),237:238)] %>%
  tidyr::pivot_longer(cols=grep("fdr",colnames(reg), value=T), names_to="superfine.cell.class", values_to="fdr", names_prefix="fdr_")
logfc.df = reg[,c(grep("logFC",colnames(reg)),237:238)] %>%
  tidyr::pivot_longer(cols=grep("logFC",colnames(reg), value=T), names_to="superfine.cell.class", values_to="logfc", names_prefix="logFC_")

enrich.df = left_join(tstats.df, fdr.df, by=c("superfine.cell.class","ensembl","gene")) %>%
  left_join(logfc.df, by=c("superfine.cell.class","ensembl","gene"))
enrich.df = left_join(enrich.df, 
                      distinct(as.data.frame(colData(sce)[,c("fine.cell.class","superfine.cell.class")])) %>%
                        mutate(superfine.cell.class=gsub("/",".",superfine.cell.class)),
                      by="superfine.cell.class")

deep.df = filter(enrich.df, superfine.cell.class %in% c("Sub.1","Sub.2","L6.1","L6.2","L6b","L5.1","L5.2"),
                 fdr<.0001 & logfc>2)

top50.40 = names(sort(loads[,"nmf40"], decreasing=T))[1:50]
filter(deep.df, gene %in% top50.40) %>% group_by(gene) %>% add_tally() %>%
  mutate(sig.classes= paste(superfine.cell.class, collapse=" ")) %>%
  filter(superfine.cell.class=="Sub.1") %>% dplyr::select(gene, sig.classes, t_stat, fdr, logfc, n) %>%
  arrange(n)
#"PEX5L","ENOX1","CDH8","GRIK1","TOX"

top50.54 = names(sort(loads[,"nmf54"], decreasing=T))[1:50]
filter(deep.df, gene %in% top50.54) %>% group_by(gene) %>% add_tally() %>%
  mutate(sig.classes= paste(superfine.cell.class, collapse=" ")) %>%
  filter(superfine.cell.class=="Sub.2") %>% dplyr::select(gene, sig.classes, t_stat, fdr, logfc, n) %>%
  arrange(n)
#"ZNF385D","CNTN4","LUZP2","TSHZ2","SEMA5A","TOX"

top50.65 = names(sort(loads[,"nmf65"], decreasing=T))[1:50]
filter(deep.df, gene %in% top50.65) %>% group_by(gene) %>% add_tally() %>%
  mutate(sig.classes= paste(superfine.cell.class, collapse=" ")) %>%
  filter(superfine.cell.class %in% c("L6.1")) %>% dplyr::select(gene, sig.classes, t_stat, fdr, logfc, n) %>%
  arrange(n)
#just 1: AL391117.1

top50.53 = names(sort(loads[,"nmf53"], decreasing=T))[1:50]
filter(deep.df, gene %in% top50.53) %>% group_by(gene) %>% add_tally() %>%
  mutate(sig.classes= paste(superfine.cell.class, collapse=" ")) %>%
  filter(superfine.cell.class %in% c("L6b")) %>% dplyr::select(gene, sig.classes, t_stat, fdr, logfc, n) %>%
  arrange(n)
#"SGCZ","MARCH1","LRMDA","NFIA","TENM3"

#dotplot for supp
plotDots(sce_deep, features=c("CDH8","GRIK1","PEX5L","ENOX1","TOX",#nmf40
                              "CNTN4","TSHZ2","ZNF385D","SEMA5A","LUZP2",#nmf54
                              "AL391117.1",#nmf65
                              "TENM3","SGCZ","NFIA","MARCH1","LRMDA"#nmf53
                              ), group="new.classes")+
  scale_color_viridis_c(option="F", direction=-1)+
  labs(title="Top nmf genes that are DEGs that support L6.1 as deep sub")+
  theme(axis.text.x=element_text(angle=45, hjust=1), text=element_text(size=16))


##############################
######### cortical layer markers plots - main figure
##############################
sce_subset$plot.cell.class = as.character(sce_subset$superfine.cell.class)
sce_subset$plot.cell.class = ifelse(sce_subset$fine.cell.class %in% c("HATA","Amy","Thal","Cajal","GC","CA2-4"), 
                                    as.character(sce_subset$fine.cell.class), sce_subset$plot.cell.class)
sce_subset$new.cell.class = factor(sce_subset$plot.cell.class,
                                   levels=c("GC","CA2-4","HATA","Amy","Thal","Cajal","L2/3.5",
                                            "CA1","ProS",
                                            "Sub.1","Sub.2","L6.1","L2/3.1",
                                            "L6b","L6.2","L5.2","L5.1",
                                            "L2/3.2","L2/3.3","L2/3.6","L2/3.4"),
                                   labels=c("GC","CA2-4","HATA","Amy","Thal","Cajal","L2/3.5",
                                            "CA1","ProS",
                                            "Sub.1","Sub.2","Sub.3","PaS",
                                            "PAC.L6b","PAC.L6","PAC.CBLN2+","ENT.L5",
                                            "ENT.sup4","ENT.sup3","ENT.sup2","ENT.sup1")
)

new.class.palette = c("GC"="grey50","CA2-4"="grey50","HATA"="grey50","Amy"="grey50","Thal"="grey50","Cajal"="grey50","L2/3.5"="grey50",
                       "CA1"="#984ea3","ProS"="#e41a1c","Sub.1"="#ff7f00","Sub.2"="#fdbf6f","Sub.3"="#fc4e2a","PaS"="#f768a1",
                       "PAC.L6b"="#023858","PAC.L6"="#045a8d","ENT.L5"="#016c59","PAC.CBLN2+"="#a6bddb",
                       "ENT.sup4"="#006837","ENT.sup3"="#78c679","ENT.sup2"="#00bb00","ENT.sup1"="#add294")
#classic cortical markers
plotExpression(sce_subset[,!sce_subset$new.cell.class %in% c("Thal","Cajal")], 
	       features=c("SATB2","TLE4","CUX2","RORB","BCL11B","CBLN2"), x="new.cell.class", 
               color_by = "new.cell.class", point_size=.5)+
  scale_color_manual(values=new.class.palette)+
  facet_wrap(vars(Feature), ncol=2)+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), 
                                          legend.position="none",text=element_text(size=14),
                                          axis.title.x=element_blank())
#DEGs from Fig3C
plotExpression(sce_subset[,!sce_subset$new.cell.class %in% c("Thal","Cajal")], 
               features=c("FN1","COL24A1","POU3F1","PART1","KCNH5","TESPA1"), x="new.cell.class", 
               color_by = "new.cell.class", point_size=.5)+
  scale_color_manual(values=new.class.palette)+
  facet_wrap(vars(Feature), ncol=2, scales="free_y")+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), 
                                          legend.position="none",text=element_text(size=14),
                                          axis.title.x=element_blank())

##############################
######### subiculum-focused DE
##############################

#subiculum DE
sce_sub = sce_subset[,sce_subset$superfine.cell.class %in% c("CA1","Sub.1","Sub.2","ProS","L6.1","L2/3.1")]
sce_sub$superfine.cell.class = droplevels(sce_sub$superfine.cell.class)
table(rowSums(counts(sce_sub)>0)>100)
sce_sub = sce_sub[rowSums(counts(sce_sub)>0)>100,]
sub.results = scran::findMarkers(sce_sub, group=sce_sub$superfine.cell.class, test="binom")

filter(as.data.frame(sub.results[['Sub.1']]), -log10(FDR)>30 & logFC.CA1>1) %>% nrow() #897
filter(as.data.frame(sub.results[['Sub.1']]), -log10(FDR)>30 & logFC.ProS>1) %>% nrow() #666
filter(as.data.frame(sub.results[['Sub.1']]), -log10(FDR)>30 & logFC.Sub.2>1) %>% nrow() #670
filter(as.data.frame(sub.results[['Sub.1']]), -log10(FDR)>30 & logFC.L6.1>1) %>% nrow() #1198
filter(as.data.frame(sub.results[['Sub.1']]), -log10(FDR)>30 & logFC.L2.3.1>1) %>% nrow() #1730

sub.1.list <- list("CA1"=rownames(filter(as.data.frame(sub.results[['Sub.1']]), -log10(FDR)>30 & logFC.CA1>1)),
                   "ProS"=rownames(filter(as.data.frame(sub.results[['Sub.1']]), -log10(FDR)>30 & logFC.ProS>1)),
                   "Sub.2"=rownames(filter(as.data.frame(sub.results[['Sub.1']]), -log10(FDR)>30 & logFC.Sub.2>1)),
                   "L6.1"=rownames(filter(as.data.frame(sub.results[['Sub.1']]), -log10(FDR)>30 & logFC.L6.1>1)),
                   "L2/3.1"=rownames(filter(as.data.frame(sub.results[['Sub.1']]), -log10(FDR)>30 & logFC.L2.3.1>1)))
sub.2.list <- list("CA1"=rownames(filter(as.data.frame(sub.results[['Sub.2']]), -log10(FDR)>30 & logFC.CA1>1)),
                   "ProS"=rownames(filter(as.data.frame(sub.results[['Sub.2']]), -log10(FDR)>30 & logFC.ProS>1)),
                   "Sub.1"=rownames(filter(as.data.frame(sub.results[['Sub.2']]), -log10(FDR)>30 & logFC.Sub.1>1)),
                   "L6.1"=rownames(filter(as.data.frame(sub.results[['Sub.2']]), -log10(FDR)>30 & logFC.L6.1>1)),
                   "L2/3.1"=rownames(filter(as.data.frame(sub.results[['Sub.2']]), -log10(FDR)>30 & logFC.L2.3.1>1)))
pros.list <- list("CA1"=rownames(filter(as.data.frame(sub.results[['ProS']]), -log10(FDR)>30 & logFC.CA1>1)),
                   "Sub.2"=rownames(filter(as.data.frame(sub.results[['ProS']]), -log10(FDR)>30 & logFC.Sub.2>1)),
                   "Sub.1"=rownames(filter(as.data.frame(sub.results[['ProS']]), -log10(FDR)>30 & logFC.Sub.1>1)),
                   "L6.1"=rownames(filter(as.data.frame(sub.results[['ProS']]), -log10(FDR)>30 & logFC.L6.1>1)),
                   "L2/3.1"=rownames(filter(as.data.frame(sub.results[['ProS']]), -log10(FDR)>30 & logFC.L2.3.1>1)))
#sub.1
t1 = table(unlist(sub.1.list))
sub1.unique = names(t1)[t1==length(sub.1.list)]
length(sub1.unique) #146
avg.expr = mutate(as.data.frame(sub.results[['Sub.1']][sub1.unique,]), avg.sub1.expr = rowMeans(logcounts(sce_sub)[sub1.unique,sce_sub$superfine.cell.class=="Sub.1"]))
rownames(filter(avg.expr, avg.sub1.expr>1))
#[1] "AC007368.1" "AL138694.1" "ATP6V1C2"   "COL21A1"    "EBF4"      
#[6] "FN1"        "NDST4"      "PARD3B"     "PRKCH"      "RAPGEF3"
plotExpression(sce_sub, features = c("COL21A1","NDST4","EBF4","ATP6V1C2","PRKCH","RAPGEF3"), 
               x="new.cell.class", exprs_values = "logcounts", colour_by = "new.cell.class")+
  scale_color_manual(values=new.class.palette)+theme_bw()+
  theme(text=element_text(size=16), legend.position="none", axis.title.x=element_blank(), 
        axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
#sub.2
t2 = table(unlist(sub.2.list))
sub2.unique = names(t2)[t2==length(sub.2.list)]
length(sub2.unique) #184
avg.expr = mutate(as.data.frame(sub.results[['Sub.2']][sub2.unique,]), 
                  avg.sub2.expr = rowMeans(logcounts(sce_sub)[sub2.unique,sce_sub$superfine.cell.class=="Sub.2"]))
rownames(filter(avg.expr, avg.sub2.expr>1))
#"GDNF-AS1" "LHFPL3"   "MAMDC2"   "PCED1B"   "RORB"     "SULF1"    "TRPC3" 
plotExpression(sce_sub, features = c("GDNF-AS1","LHFPL3","MAMDC2","PCED1B","SULF1","TRPC3"), 
               x="new.cell.class", exprs_values = "logcounts", colour_by = "new.cell.class")+
  scale_color_manual(values=new.class.palette)+theme_bw()+
  theme(text=element_text(size=16), legend.position="none", axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
#pros
t3 = table(unlist(pros.list))
pros.unique = names(t3)[t3==length(pros.list)]
length(pros.unique) #17
avg.expr = mutate(as.data.frame(sub.results[['ProS']][pros.unique,]), 
                  avg.pros.expr = rowMeans(logcounts(sce_sub)[pros.unique,sce_sub$superfine.cell.class=="ProS"]))
rownames(filter(avg.expr, avg.pros.expr>1)) 
#none

#deep DE
sce_subd = sce_subset[,sce_subset$superfine.cell.class %in% c("Sub.1","Sub.2","L6.1","L6b","L6.2")]
sce_subd$superfine.cell.class = droplevels(sce_subd$superfine.cell.class)
sce_subd = sce_subd[rowSums(counts(sce_subd)>0)>100,]
subd.results = scran::findMarkers(sce_subd, group=sce_subd$superfine.cell.class, test="binom")

filter(as.data.frame(subd.results[['L6.1']]), -log10(FDR)>30 & logFC.Sub.1>1) %>% nrow() #198
filter(as.data.frame(subd.results[['L6.1']]), -log10(FDR)>30 & logFC.Sub.2>1) %>% nrow() #119
filter(as.data.frame(subd.results[['L6.1']]), -log10(FDR)>30 & logFC.L6b>1) %>% nrow() #155
filter(as.data.frame(subd.results[['L6.1']]), -log10(FDR)>30 & logFC.L6.2>1) %>% nrow() #184

l61.list <- list("L6b"=rownames(filter(as.data.frame(subd.results[['L6.1']]), -log10(FDR)>30 & logFC.L6b>1)),
                  "Sub.2"=rownames(filter(as.data.frame(subd.results[['L6.1']]), -log10(FDR)>30 & logFC.Sub.2>1)),
                  "Sub.1"=rownames(filter(as.data.frame(subd.results[['L6.1']]), -log10(FDR)>30 & logFC.Sub.1>1)),
                  "L6.2"=rownames(filter(as.data.frame(subd.results[['L6.1']]), -log10(FDR)>30 & logFC.L6.2>1)))

t4 = table(unlist(l61.list))
l61.unique = names(t4)[t4==length(l61.list)]
length(l61.unique) #68
avg.expr = mutate(as.data.frame(subd.results[['L6.1']][l61.unique,]), 
                  avg.l61.expr = rowMeans(logcounts(sce_subd)[l61.unique,sce_subd$superfine.cell.class=="L6.1"]))
rownames(filter(avg.expr, avg.l61.expr>1)) #
#[1] "AC007100.1" "AC010967.1" "AC023503.1" "AC046195.2" "AL356295.1"
#[6] "CD36"       "COL4A1"     "COL4A2"     "DISC1"      "FSTL5"     
#[11] "GUCA1C"     "LINC01194"  "LINC01239"  "LINC01821"  "NR2F2-AS1" 
#[16] "PLEKHG1"    "RASGEF1B"   "SCN7A"      "SCUBE1"     "SNCAIP"    
#[21] "VEGFC" 
plotExpression(sce_subd, features = c(#"COL4A1","COL4A2","VEGFC",
  "FSTL5","PLEKHG1","RASGEF1B","GUCA1C","SCN7A","SCUBE1","SNCAIP","DISC1"), 
  x="new.cell.class", exprs_values = "logcounts", colour_by = "new.cell.class")+
  scale_color_manual(values=new.class.palette)+theme_bw()+
  theme(text=element_text(size=16), legend.position="none", axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

#superficial DE
sce_subs = sce_subset[,sce_subset$superfine.cell.class %in% c("Sub.1","Sub.2","ProS","L2/3.1","L2/3.2","L2/3.4","L2/3.3","L2/3.6")]
sce_subs$superfine.cell.class = droplevels(sce_subs$superfine.cell.class)
sce_subs = sce_subs[rowSums(counts(sce_subs)>0)>100,]
subs.results = scran::findMarkers(sce_subs, group=sce_subs$superfine.cell.class, test="binom")

filter(as.data.frame(subs.results[['L2/3.1']]), -log10(FDR)>30 & logFC.Sub.1>1) %>% nrow() #446
filter(as.data.frame(subs.results[['L2/3.1']]), -log10(FDR)>30 & logFC.Sub.2>1) %>% nrow() #419
filter(as.data.frame(subs.results[['L2/3.1']]), -log10(FDR)>30 & logFC.ProS>1) %>% nrow() #463
filter(as.data.frame(subs.results[['L2/3.1']]), -log10(FDR)>30 & logFC.L2.3.2>1) %>% nrow() #351
filter(as.data.frame(subs.results[['L2/3.1']]), -log10(FDR)>30 & logFC.L2.3.3>1) %>% nrow() #418
filter(as.data.frame(subs.results[['L2/3.1']]), -log10(FDR)>30 & logFC.L2.3.4>1) %>% nrow() #321
filter(as.data.frame(subs.results[['L2/3.1']]), -log10(FDR)>30 & logFC.L2.3.6>1) %>% nrow() #362

l23.list <- list("ProS"=rownames(filter(as.data.frame(subs.results[['L2/3.1']]), -log10(FDR)>30 & logFC.ProS>1)),
                 "Sub.2"=rownames(filter(as.data.frame(subs.results[['L2/3.1']]), -log10(FDR)>30 & logFC.Sub.2>1)),
                 "Sub.1"=rownames(filter(as.data.frame(subs.results[['L2/3.1']]), -log10(FDR)>30 & logFC.Sub.1>1)),
                 "L2/3.2"=rownames(filter(as.data.frame(subs.results[['L2/3.1']]), -log10(FDR)>30 & logFC.L2.3.2>1)),
                 "L2/3.3"=rownames(filter(as.data.frame(subs.results[['L2/3.1']]), -log10(FDR)>30 & logFC.L2.3.3>1)),
                 "L2/3.4"=rownames(filter(as.data.frame(subs.results[['L2/3.1']]), -log10(FDR)>30 & logFC.L2.3.4>1)),
                 "L2/3.6"=rownames(filter(as.data.frame(subs.results[['L2/3.1']]), -log10(FDR)>30 & logFC.L2.3.6>1)))

t5 = table(unlist(l23.list))
l23.unique = names(t5)[t5==length(l23.list)]
length(l23.unique) #38
avg.expr = mutate(as.data.frame(subs.results[['L2/3.1']][l23.unique,]), 
                  avg.l23.expr = rowMeans(logcounts(sce_subs)[l23.unique,sce_subs$superfine.cell.class=="L2/3.1"]))
rownames(filter(avg.expr, avg.l23.expr>1)) #5
#[1] "AC008662.1" "AL161629.1" "FSTL5"      "MDFIC"      "WSCD1" 
plotExpression(sce_subs, features = c("AC008662.1","AL161629.1","FSTL5","MDFIC","WSCD1"),
               x="new.cell.class", exprs_values = "logcounts", colour_by = "new.cell.class")+
  scale_color_manual(values=new.class.palette)+theme_bw()+
  theme(text=element_text(size=16), legend.position="none", axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

#combine all the new DEGs into a dot plot - main figure
plotDots(sce_subset[,!sce_subset$new.cell.class %in% c("Thal","Cajal")], 
         features=c("COL21A1","NDST4","EBF4","ATP6V1C2","PRKCH","RAPGEF3",
                    "GDNF-AS1","LHFPL3","MAMDC2","PCED1B","TRPC3",#"RORB",
                    "PLEKHG1","RASGEF1B","GUCA1C","SCN7A","SCUBE1",#"SNCAIP",
                    "FSTL5","AC008662.1","AL161629.1","MDFIC","WSCD1"),
         group = "new.cell.class")+
  scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        text=element_text(size=14), axis.title.x=element_blank())
