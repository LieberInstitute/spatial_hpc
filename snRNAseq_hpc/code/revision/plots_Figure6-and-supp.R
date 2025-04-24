library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scater)

set.seed(123)

load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))

## tsne
sce_pyr = sce[,sce$fine.cell.class %in% c("CA1/ProS","L2/3.PrS.PaS","L2/3.PrS.Ent","L5/6","L6/6b","Sub.1","Sub.2")]
reducedDim(sce_pyr, "TSNE") = calculateTSNE(sce_pyr, dimred="MNN", n_dimred=50)

nmf.list = c("nmf15","nmf40","nmf54","nmf84","nmf45","nmf27","nmf51","nmf68","nmf22","nmf53","nmf65","nmf17")

plist = lapply(nmf.list, function(x) 
  plotReducedDim(sce_pyr, dimred="TSNE", color_by=x, point_size=.3)+
    scale_color_viridis_c(option="F", direction=-1)+ggtitle(x)+
    theme(legend.position="none", axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(),
          plot.title=element_text(hjust=.5, size=14))
)

ggsave("snRNAseq_hpc/plots/revision/Figure6_nmf-tsne.png",do.call(gridExtra::grid.arrange, c(plist, ncol=4)),
       bg="white", height=6, width=7, units="in")

##############################
######### L6.1 as part of subiculum
##############################

sce_subset = sce[,sce$mid.cell.class=="ExcN"]
sce_subset$superfine.cell.class = droplevels(sce_subset$superfine.cell.class)

tmp = sce_subset[,sce_subset$superfine.cell.class=="L6b"]
p0 <- ggplot(as.data.frame(colData(sce_tmp)),
       aes(x=nmf53, y=nmf65, color=nmf65>.00035))+
  geom_point()+theme_bw()+ggtitle("L6b only")+
  scale_color_manual(values=c("black","red3"))
ggsave(file="snRNAseq_hpc/plots/revision/ED_L6b-nmf65_scatter.pdf", p0, height=6, width=6)
tmp$nmf65pos = tmp$nmf65>.00035


#sce_subset$is_deep = ifelse(sce_subset$fine.cell.class %in% c("L5/6","L6/6b","Sub.1","Sub.2"), sce_subset$superfine.cell.class, "other")
sce_deep = sce_subset[,sce_subset$fine.cell.class %in% c("L5/6","L6/6b","Sub.1","Sub.2")]
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

try.65 = setdiff(top50.65, c(top50.40, top50.54, top50.53))
#[1] "MGAT4C"     "CTNND2"     "SLC35F3"    "CDH13"      "KHDRBS3"   
#[6] "RIMS2"      "NAV3"       "AL391117.1" "GALNT17"    "NEGR1"     
#[11] "MMP16"      "PID1"       "UTRN"       "TMEM108"    "ADGRL3" 
filter(deep.df, gene %in% try.65) %>% group_by(gene) %>% add_tally() %>%
  mutate(sig.classes= paste(superfine.cell.class, collapse=" ")) %>%
  filter()

plotDots(sce_deep, features=try.65, group="new.classes")+
  scale_color_gradient(low="grey90", high="black")+
  labs(size="prop. nuclei", color="Avg. expr.")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), text=element_text(size=11))
#add UTRN


#dotplot for supp
p1 <- plotDots(sce_deep, features=c("CDH8","GRIK1","PEX5L","ENOX1","TOX",#nmf40
                              "CNTN4","TSHZ2","ZNF385D","SEMA5A","LUZP2",#nmf54
                              "AL391117.1","UTRN",#nmf65
                              "TENM3","SGCZ","NFIA","MARCH1","LRMDA"#nmf53
), group="new.classes")+
  scale_color_gradient(low="grey90", high="black")+
  labs(size="prop. nuclei", color="Avg. expr.")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), text=element_text(size=11))

ggsave(file="snRNAseq_hpc/plots/revision/ED_L6b-nmf65_dotplot.pdf", p1,
       height=8, width=6, units="in")

#additional tsne for supp with clusters labeled
deep.palette = c("Sub.1"="#ff7f00","Sub.2"="#fdbf6f","L6.1"="#fc4e2a",
                      "L6b"="#023858","L6.2"="#045a8d","L5.1"="#016c59","L5.2"="#a6bddb","other"="grey")
sce_pyr$deep_clus = ifelse(sce_pyr$superfine.cell.class %in% c("Sub.1","Sub.2","L6.1","L6.2","L6b","L5.1","L5.2"), 
                           as.character(sce_pyr$superfine.cell.class), "other")

p2 <- plotReducedDim(sce_pyr, dimred="TSNE", color_by="deep_clus", point_size=.3)+# by.assay.type = "logcounts")+
    scale_color_manual(values=deep.palette)+
    theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank())

ggsave(file="snRNAseq_hpc/plots/revision/ED_deep-sub_clusters_tsne.pdf",
       ggrastr::rasterize(p2,layer='point',dpi=350),
       height=8, width=6, units="in")

#tsne for supp
gene.list <- c("GRIK1","TOX","ZNF385D","UTRN","LRMDA","NFIA")
plist = lapply(gene.list, function(x) {
  tmp <- plotReducedDim(sce_pyr, dimred="TSNE", color_by=x, point_size=.3, by.assay.type = "logcounts")+
    scale_color_gradient(low="grey90", high="black")+labs(color="")+ggtitle(x)+
    theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(),
          plot.title=element_text(hjust=.5, size=11))
  ggrastr::rasterize(tmp,layer='point',dpi=350)
})
do.call(gridExtra::grid.arrange, c(plist, ncol=2))
ggsave(file="snRNAseq_hpc/plots/revision/ED_deep-sub_markers_tsne.pdf",
       do.call(gridExtra::grid.arrange, c(plist, ncol=2)),
       height=8, width=6, units="in")

#another tsne for supp
nmf.list <- c("nmf40","nmf54","nmf68","nmf22","nmf53","nmf65")
plist2 = lapply(nmf.list, function(x) {
  tmp <- plotReducedDim(sce_pyr, dimred="TSNE", color_by=x, point_size=.3)+
    scale_color_viridis_c(option="F", direction=-1)+labs(color="")+ggtitle(x)+
    theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(),
          plot.title=element_text(hjust=.5, size=11))
  ggrastr::rasterize(tmp,layer='point',dpi=350)
}
)
ggsave(file="snRNAseq_hpc/plots/revision/ED_deep-sub_nmfs_tsne.pdf",
       do.call(gridExtra::grid.arrange, c(plist2, ncol=2)),
       height=8, width=6, units="in")


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
                                            "Sub.1","Sub.2","Sub.3","PreS",
                                            "RHP.L6b","RHP.L6","RHP.CBLN2+","ENT.L5",
                                            "ENT.sup3","ENT.sup2b","ENT.sup2a","ENT.sup1")
)

new.class.palette = c("GC"="grey50","CA2-4"="grey50","HATA"="grey50","Amy"="grey50","Thal"="grey50","Cajal"="grey50","L2/3.5"="grey50",
                      "CA1"="#984ea3","ProS"="#e41a1c","Sub.1"="#ff7f00","Sub.2"="#fdbf6f","Sub.3"="#fc4e2a","PreS"="#f768a1",
                      "RHP.L6b"="#023858","RHP.L6"="#045a8d","ENT.L5"="#016c59","RHP.CBLN2+"="#a6bddb",
                      "ENT.sup3"="#006837","ENT.sup2b"="#78c679","ENT.sup2a"="#78c679","ENT.sup1"="#add294")

#classic cortical markers and new sub markers
ctx.markers = c("SATB2","TLE4","CUX2","FN1","COL24A1","TOX")

vln.df = cbind.data.frame("new.cell.class"=colData(sce_subset)[,c("new.cell.class")], as.matrix(t(logcounts(sce_subset)[ctx.markers,])))
long.df = filter(vln.df, !new.cell.class %in% c("Thal", "Cajal")) %>% 
  tidyr::pivot_longer(all_of(ctx.markers), names_to="marker.gene", values_to="expr") %>%
  mutate(marker.gene=factor(marker.gene, levels=ctx.markers))


vln <- ggplot(long.df, aes(x=new.cell.class, y=expr, fill=new.cell.class))+
  geom_violin(scale="width")+scale_fill_manual(values=new.class.palette)+
  facet_wrap(vars(marker.gene), ncol=2)+
  theme_minimal()+labs(y="Expression (logcounts)", x="")+
  theme(legend.position="none", text=element_text(size=12), panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=.5), axis.title.x=element_blank(),
        strip.text=element_text(size=12),
        strip.background=element_rect(fill="white", color="transparent"))

ggsave("snRNAseq_hpc/plots/revision/Figure6_violin.pdf", vln, bg="white", height=8, width=7, units="in")

##############################
######### subiculum-focused DE
##############################

plot.list = readRDS("snRNAseq_hpc/processed-data/revision/subiculum-DE_sig-genes-list.rda")
plot.list
#$Sub.1
#[1] "AC007368.1" "AL138694.1" "ATP6V1C2"   "COL21A1"    "EBF4"      
#[6] "FN1"        "NDST4"      "PARD3B"     "PRKCH"      "RAPGEF3"   

#$Sub.2
#[1] "GDNF-AS1" "LHFPL3"   "MAMDC2"   "PCED1B"   "RORB"     "SULF1"    "TRPC3"   

#$ProS
#character(0)

#$Sub.3
#[1] "AC007100.1" "AC010967.1" "AC023503.1" "AC046195.2" "AL356295.1"
#[6] "CD36"       "COL4A1"     "COL4A2"     "DISC1"      "FSTL5"     
#[11] "GUCA1C"     "LINC01194"  "LINC01239"  "LINC01821"  "NR2F2-AS1" 
#[16] "PLEKHG1"    "RASGEF1B"   "SCN7A"      "SCUBE1"     "SNCAIP"    
#[21] "VEGFC"     

#$PreS
#[1] "AC008662.1" "AL161629.1" "FSTL5"      "MDFIC"      "WSCD1" 

#subiculum DE
sce_sub = sce_subset[,sce_subset$new.cell.class %in% c("CA1","Sub.1","Sub.2","ProS","Sub.3","PreS")]
sce_sub$new.cell.class = droplevels(sce_sub$new.cell.class)


plotExpression(sce_sub, features = c("COL21A1","NDST4","EBF4","ATP6V1C2","PRKCH","RAPGEF3"), 
               x="new.cell.class", exprs_values = "logcounts", colour_by = "new.cell.class")+
  scale_color_manual(values=new.class.palette)+theme_bw()+
  theme(text=element_text(size=16), legend.position="none", axis.title.x=element_blank(), 
        axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
 
m1 = c("COL21A1","NDST4","EBF4","ATP6V1C2","PRKCH","RAPGEF3")
vln.df = cbind.data.frame("new.cell.class"=colData(sce_sub)[,c("new.cell.class")], as.matrix(t(logcounts(sce_sub)[m1,])))
long.df = tidyr::pivot_longer(vln.df, all_of(m1), names_to="marker.gene", values_to="expr") %>%
  mutate(marker.gene=factor(marker.gene, levels=m1))


vln <- ggplot(long.df, aes(x=new.cell.class, y=expr, fill=new.cell.class))+
  geom_violin(scale="width")+scale_fill_manual(values=new.class.palette)+
  facet_wrap(vars(marker.gene), ncol=2)+
  theme_minimal()+labs(y="Expression (logcounts)", x="")+
  theme(legend.position="none", text=element_text(size=12), panel.grid.minor=element_blank(),
        panel.border = element_rect(fill="transparent", color="black", linewidth=.5),
        axis.text.x=element_text(angle=90, hjust=1, vjust=.3), axis.title.x=element_blank(),
        strip.text=element_text(size=12),
        strip.background=element_rect(fill="white", color="transparent"))

ggsave("snRNAseq_hpc/plots/revision/ED_Sub.1-markers_violin.pdf", vln, bg="white", height=8, width=7, units="in")


plotExpression(sce_sub, features = c("GDNF-AS1","LHFPL3","MAMDC2","PCED1B","SULF1","TRPC3"), 
               x="new.cell.class", exprs_values = "logcounts", colour_by = "new.cell.class")+
  scale_color_manual(values=new.class.palette)+theme_bw()+
  theme(text=element_text(size=16), legend.position="none", axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

m1 = c("GDNF-AS1","LHFPL3","MAMDC2","PCED1B","SULF1","TRPC3")
vln.df = cbind.data.frame("new.cell.class"=colData(sce_sub)[,c("new.cell.class")], as.matrix(t(logcounts(sce_sub)[m1,])))
long.df = tidyr::pivot_longer(vln.df, all_of(m1), names_to="marker.gene", values_to="expr") %>%
  mutate(marker.gene=factor(marker.gene, levels=m1))


vln <- ggplot(long.df, aes(x=new.cell.class, y=expr, fill=new.cell.class))+
  geom_violin(scale="width")+scale_fill_manual(values=new.class.palette)+
  facet_wrap(vars(marker.gene), ncol=2)+
  theme_minimal()+labs(y="Expression (logcounts)", x="")+
  theme(legend.position="none", text=element_text(size=12), panel.grid.minor=element_blank(),
        panel.border = element_rect(fill="transparent", color="black", linewidth=.5),
        axis.text.x=element_text(angle=90, hjust=1, vjust=.3), axis.title.x=element_blank(),
        strip.text=element_text(size=12),
        strip.background=element_rect(fill="white", color="transparent"))

ggsave("snRNAseq_hpc/plots/revision/ED_Sub.2-markers_violin.pdf", vln, bg="white", height=8, width=7, units="in")


#deep DE
sce_subd = sce_subset[,sce_subset$new.cell.class %in% c("Sub.1","Sub.2","Sub.3","RHP.L6b","RHP.L6")]
sce_subd$new.cell.class = droplevels(sce_subd$new.cell.class)

plotExpression(sce_subd, features = c(#"COL4A1","COL4A2","VEGFC",
  "FSTL5","PLEKHG1","RASGEF1B","GUCA1C","SCN7A","SCUBE1","SNCAIP","DISC1"), 
  x="new.cell.class", exprs_values = "logcounts", colour_by = "new.cell.class")+
  scale_color_manual(values=new.class.palette)+theme_bw()+
  theme(text=element_text(size=16), legend.position="none", axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

m1 = c("FSTL5","PLEKHG1","RASGEF1B","GUCA1C","SCN7A","SCUBE1","SNCAIP","DISC1")
vln.df = cbind.data.frame("new.cell.class"=colData(sce_subd)[,c("new.cell.class")], as.matrix(t(logcounts(sce_subd)[m1,])))
long.df = tidyr::pivot_longer(vln.df, all_of(m1), names_to="marker.gene", values_to="expr") %>%
  mutate(marker.gene=factor(marker.gene, levels=m1))


vln <- ggplot(long.df, aes(x=new.cell.class, y=expr, fill=new.cell.class))+
  geom_violin(scale="width")+scale_fill_manual(values=new.class.palette)+
  facet_wrap(vars(marker.gene), ncol=2)+
  theme_minimal()+labs(y="Expression (logcounts)", x="")+
  theme(legend.position="none", text=element_text(size=12), panel.grid.minor=element_blank(),
        panel.border = element_rect(fill="transparent", color="black", linewidth=.5),
        axis.text.x=element_text(angle=90, hjust=1, vjust=.3), axis.title.x=element_blank(),
        strip.text=element_text(size=12),
        strip.background=element_rect(fill="white", color="transparent"))

ggsave("snRNAseq_hpc/plots/revision/ED_Sub.3-markers_violin.pdf", vln, bg="white", height=8, width=7, units="in")



#superficial DE
sce_subs = sce_subset[,sce_subset$new.cell.class %in% c("Sub.1","Sub.2","ProS","PreS","ENT.sup1","ENT.sup2a","ENT.sup2b","ENT.sup3")]
sce_subs$new.cell.class = droplevels(sce_subs$new.cell.class)

plotExpression(sce_subs, features = c("AC008662.1","AL161629.1","FSTL5","MDFIC","WSCD1"),
               x="new.cell.class", exprs_values = "logcounts", colour_by = "new.cell.class")+
  scale_color_manual(values=new.class.palette)+theme_bw()+
  theme(text=element_text(size=16), legend.position="none", axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

m1 = c("AC008662.1","AL161629.1","FSTL5","MDFIC","WSCD1")
vln.df = cbind.data.frame("new.cell.class"=colData(sce_subs)[,c("new.cell.class")], as.matrix(t(logcounts(sce_subs)[m1,])))
long.df = tidyr::pivot_longer(vln.df, all_of(m1), names_to="marker.gene", values_to="expr") %>%
  mutate(marker.gene=factor(marker.gene, levels=m1))


vln <- ggplot(long.df, aes(x=new.cell.class, y=expr, fill=new.cell.class))+
  geom_violin(scale="width")+scale_fill_manual(values=new.class.palette)+
  facet_wrap(vars(marker.gene), ncol=2)+
  theme_minimal()+labs(y="Expression (logcounts)", x="")+
  theme(legend.position="none", text=element_text(size=12), panel.grid.minor=element_blank(),
        panel.border = element_rect(fill="transparent", color="black", linewidth=.5),
        axis.text.x=element_text(angle=90, hjust=1, vjust=.3), axis.title.x=element_blank(),
        strip.text=element_text(size=12),
        strip.background=element_rect(fill="white", color="transparent"))

ggsave("snRNAseq_hpc/plots/revision/ED_PreS-markers_violin.pdf", vln, bg="white", height=8, width=7, units="in")



#combine all the new DEGs into a dot plot - main figure
dplot <- plotDots(sce_subset[,!sce_subset$new.cell.class %in% c("Thal","Cajal")], 
         features=c("COL21A1","NDST4","EBF4","ATP6V1C2","PRKCH","RAPGEF3",
                    "GDNF-AS1","LHFPL3","MAMDC2","PCED1B","TRPC3",#"RORB",
                    "PLEKHG1","RASGEF1B","GUCA1C","SCN7A","SCUBE1",#"SNCAIP",
                    "FSTL5","AC008662.1","AL161629.1","MDFIC","WSCD1"),
         group = "new.cell.class")+
  scale_color_gradient(low="white", high="black")+
  scale_size(limits=c(0,1))+
  labs(size="Prop. nuclei",color="Avg. expr.")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
	axis.text.y=element_text(size=12), axis.title.y=element_blank(),
        text=element_text(size=14), axis.title.x=element_blank())

ggsave("snRNAseq_hpc/plots/revision/Figure6_dotplot.pdf", dplot, bg="white", height=8, width=7, units="in")
