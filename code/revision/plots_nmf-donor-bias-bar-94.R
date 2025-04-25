library(SingleCellExperiment)
library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(scater)
library(ggspavis)

set.seed(123)
#load datasets
load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))
load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))
#subset patterns
non0.spots = colSums(as.matrix(colData(spe)[,paste0("nmf",1:100)])>0)
remove.nmf = names(non0.spots[non0.spots<1050])
remove.nmf[c(2,3,6)] <- c("nmf2","nmf3","nmf16")
remove.nmf = c(remove.nmf, "nmf37", "nmf28")

nmf.ordered = c("nmf37", "nmf9", "nmf3", "nmf1", "nmf55", "nmf56", "nmf57", "nmf58", "nmf59", "nmf75", "nmf80", "nmf91", "nmf92", "nmf94", "nmf24",#broad/non-specific
                "nmf2", "nmf6", "nmf8", "nmf18", "nmf20", "nmf71", "nmf72", #nrn not glia
                "nmf12", "nmf28", 'nmf7',#inhb neurons
                "nmf31","nmf4", "nmf13", "nmf16", "nmf21", "nmf25", "nmf34", "nmf49", #exc nrn
                "nmf26", "nmf10", "nmf14", "nmf5", "nmf66", #GCs
                "nmf52", "nmf11", "nmf63", "nmf61", #MC, CA3, CA2
                "nmf15", "nmf32", #CA1, ProS,
                "nmf40", "nmf54", #Sub.1, Sub.2
                "nmf65", "nmf22", "nmf53", "nmf68", "nmf51", #L6, L5
                "nmf45", "nmf84", "nmf27", "nmf17", "nmf78", #L2/3
                "nmf62", "nmf69", "nmf29", "nmf43", "nmf64", #amy
                "nmf23", #thal, cajal
                "nmf67", "nmf73", "nmf83", "nmf50", "nmf35", "nmf47", "nmf88", "nmf46", "nmf60", "nmf74", "nmf86", "nmf93", #GABA 
                "nmf81", "nmf19", "nmf76", "nmf79", #astro
                "nmf42", "nmf44", "nmf38", "nmf77", "nmf33", "nmf36", #oligo
                "nmf90", "nmf39", "nmf98", "nmf82", "nmf96", "nmf100", #micro immune
                "nmf87", "nmf30", "nmf41", "nmf48", #ependy, CP
                "nmf97", "nmf89", "nmf70", "nmf85", "nmf99", "nmf95" #endo
)
nmf.ordered.keep = setdiff(nmf.ordered, remove.nmf)
nmf.ordered.remove = intersect(nmf.ordered, remove.nmf)
#nmf.ordered.remove = nmf.ordered.remove[c(1,8,2:7,9:19)] #move sex to front for 400 filter
nmf.ordered.remove = nmf.ordered.remove[c(1,10,2:9,11:34)] #move sex to front for 1050 filter

# general patterns only:
gen.nmf = c("nmf37", "nmf9", "nmf3", "nmf1", "nmf55", "nmf56", "nmf57", "nmf58", "nmf59", "nmf75", "nmf80", "nmf91", "nmf92", "nmf94", "nmf24",#broad/non-specific
            "nmf2", "nmf6", "nmf8", "nmf18", "nmf20", "nmf71", "nmf72", #nrn not glia
            "nmf12", "nmf28", 'nmf7',#inhb neurons
            "nmf31","nmf4", "nmf13", "nmf16", "nmf21", "nmf25", "nmf34", "nmf49") #exc nrn
gen.nmf <- intersect(gen.nmf, nmf.ordered.keep)

#sce
gen.non0.mtx = as.matrix(colData(sce)[,gen.nmf])
gen.non0.mtx = gen.non0.mtx>0

gen.non0 = cbind(as.data.frame(gen.non0.mtx), as.data.frame(colData(sce)[rownames(gen.non0.mtx),c("brnum","round","sort")]))
test = group_by(gen.non0, brnum, round, sort) %>% summarise_at(gen.nmf, sum) %>%
  tidyr::pivot_longer(all_of(gen.nmf), names_to="nmf", values_to="n.non0.spots")
test$nmf_f = factor(test$nmf, levels=nmf.ordered.keep)
test$xlabel = paste(test$brnum, test$sort)

sn.brnum.palette = c("#4575b4","#313695", "#b2df8a","#33a02c",
                     "#fb9a99","#e31a1c", "#ff7f00","#f46d43",
                     "#c2a5cf","#6a3d9a", "grey","grey50",
                     "#ffff99","goldenrod", "#8dd3c7","aquamarine",
                     "#f4a582","#b2182b", "grey40","black")
names(sn.brnum.palette) = paste(rep(c("Br2720","Br2743","Br3942","Br6423","Br6432",
                                        "Br6471","Br6522","Br8325","Br8492","Br8667"), each=2),
                                  rep(c("PI","PI+NeuN+"), length.out=20))

g1 <- ggplot(group_by(test, nmf_f, xlabel) %>% summarise(n.non0.spots=sum(n.non0.spots)), 
       aes(x=n.non0.spots, y=nmf_f, fill=xlabel))+
  geom_bar(stat="identity", position="fill", color="grey50", linewidth=.5, width=1)+
  scale_x_continuous(expand=c(0,0))+
  scale_fill_manual(values=sn.brnum.palette)+guides(fill = guide_legend(ncol = 2))+
  labs(x="% nuclei with non-0 loading", title="General NMF patterns: snRNA-seq",fill="Sample")+
  theme_minimal()+theme(legend.position="bottom", legend.byrow = T, axis.ticks.x = element_line(color="black"),
                        axis.title.y=element_blank(), axis.text.x=element_text(hjust=1), text=element_text(size=10))

#specific patterns only
spec.nmf = setdiff(nmf.ordered.keep, gen.nmf)
spec.non0.mtx = as.matrix(colData(sce)[,spec.nmf])
spec.non0.mtx = spec.non0.mtx>0
spec.non0 = cbind(as.data.frame(spec.non0.mtx), 
                  as.data.frame(colData(sce)[rownames(spec.non0.mtx),c("brnum","round","sort")]))
test2 = group_by(spec.non0, brnum, round, sort) %>% summarise_at(spec.nmf, sum) %>%
  tidyr::pivot_longer(all_of(spec.nmf), names_to="nmf", values_to="n.non0.spots")
test2$nmf_f = factor(test2$nmf, levels=nmf.ordered.keep)
test2$xlabel = paste(test2$brnum, test2$sort)

s1 <- ggplot(group_by(test2, nmf_f, xlabel) %>% summarise(n.non0.spots=sum(n.non0.spots)), 
       aes(x=n.non0.spots, y=nmf_f, fill=xlabel))+
  geom_bar(stat="identity", position="fill", color="grey50", linewidth=.5, width=1)+
  scale_x_continuous(expand=c(0,0))+
  scale_fill_manual(values=sn.brnum.palette)+
  labs(x="% nuclei with non-0 loading", title="Specific NMF patterns: snRNA-seq")+
  theme_minimal()+theme(legend.position="none", axis.title.y=element_blank(), axis.text.x=element_text(hjust=1),
                        axis.ticks.x = element_line(color="black"), text=element_text(size=10))

#################
## ok now do same as above but with SRT
#################

# general patterns
gen.non0.mtx = as.matrix(colData(spe)[,gen.nmf])
gen.non0.mtx = gen.non0.mtx>0
gen.non0 = cbind(as.data.frame(gen.non0.mtx), 
                 as.data.frame(colData(spe)[rownames(gen.non0.mtx),c("brnum","sample_id")]))
test3 = group_by(gen.non0, brnum, sample_id) %>% summarise_at(gen.nmf, sum) %>%
  tidyr::pivot_longer(all_of(gen.nmf), names_to="nmf", values_to="n.non0.spots")
test3$nmf_f = factor(test3$nmf, levels=nmf.ordered.keep)
test3$brnum_f = factor(test3$brnum, levels=levels(spe$brnum)[c(10,1:9)])

brnum.palette = c("#313695","#33a02c","#fb9a99","#ff7f00","#6a3d9a",
                  "grey","#ffff99","#8dd3c7","#b2182b","grey40")
names(brnum.palette) <- c("Br2720","Br2743","Br3942","Br6423","Br6432",
                          "Br6471","Br6522","Br8325","Br8492","Br8667")

g2 <- ggplot(group_by(test3, nmf_f, brnum_f) %>% summarise(n.non0.spots=sum(n.non0.spots)), 
       aes(x=n.non0.spots, y=nmf_f, fill=brnum_f))+
  geom_bar(stat="identity", position="fill", color="grey50", linewidth=.5, width=1)+
  scale_x_continuous(expand=c(0,0))+
  scale_fill_manual(values=brnum.palette)+
  guides(fill = guide_legend(ncol = 1))+
  labs(x="% spots with non-0 loading", title="General NMF patterns: SRT",fill="Donor")+
  theme_minimal()+theme(legend.position="bottom", legend.byrow = T,
                        axis.title.y=element_blank(), axis.text.x=element_text(hjust=1),
                        axis.ticks.x = element_line(color="black"), text=element_text(size=10))

#specific patterns
spec.nmf = setdiff(nmf.ordered.keep, gen.nmf)
spec.non0.mtx = as.matrix(colData(spe)[,spec.nmf])
spec.non0.mtx = spec.non0.mtx>0
spec.non0 = cbind(as.data.frame(spec.non0.mtx), 
                  as.data.frame(colData(spe)[rownames(spec.non0.mtx),c("brnum","sample_id")]))
test4 = group_by(spec.non0, brnum) %>% summarise_at(spec.nmf, sum) %>%
  tidyr::pivot_longer(all_of(spec.nmf), names_to="nmf", values_to="n.non0.spots")
test4$nmf_f = factor(test4$nmf, levels=nmf.ordered.keep)
test4$brnum_f = factor(test4$brnum, levels=levels(spe$brnum)[c(10,1:9)])

s2 <- ggplot(group_by(test4, nmf_f, brnum_f) %>% summarise(n.non0.spots=sum(n.non0.spots)), 
       aes(x=n.non0.spots, y=nmf_f, fill=brnum_f))+
  geom_bar(stat="identity", position="fill", color="grey50", linewidth=.5, width=1)+
  scale_x_continuous(expand=c(0,0))+
  scale_fill_manual(values=brnum.palette)+
  labs(x="% spots with non-0 loading", title="Specific NMF patterns: SRT")+
  theme_minimal()+theme(legend.position="none", axis.title.y=element_blank(), axis.text.x=element_text(hjust=1),
                        axis.ticks.x = element_line(color="black"), text=element_text(size=10))

pdf(file="plots/revision/supp_nmf-donor-proportion.pdf", height=12, width=8)
gridExtra::grid.arrange(g1, s1, g2, s2, ncol=2)
dev.off()  

#################
## look at nmf94 that became cell type specific
#################

#summarized dot plot
seed.94 = as.matrix(colData(sce)[,c("nmf94","nmf48")])
seed.94 = seed.94>0
d94 = cbind.data.frame(superfine.cell.class=as.data.frame(colData(sce))[,"mid.cell.class"],
                      seed.94) %>% 
  group_by(superfine.cell.class) %>% add_tally(name="total") %>%
  group_by(superfine.cell.class, total) %>%
  summarise_at(c("nmf94","nmf48"), sum) %>%
  tidyr::pivot_longer(c("nmf94","nmf48"), values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)

seed2.94 = as.matrix(colData(sce)[,c("nmf94","nmf48")])
seed2.94 = apply(seed2.94, 2, scale)
d2.94 = cbind.data.frame(superfine.cell.class=as.data.frame(colData(sce))[,"mid.cell.class"],
                      seed2.94) %>% 
  group_by(superfine.cell.class) %>%
  summarise_at(c("nmf94","nmf48"), mean) %>% 
  tidyr::pivot_longer(c("nmf94","nmf48"), values_to="scaled.avg", names_to="nmf")

sce.94.df = left_join(d94[,c("superfine.cell.class","nmf","prop")], 
                   d2.94[,c("superfine.cell.class","nmf","scaled.avg")]) %>%
  mutate(#superfine.cell.class=factor(superfine.cell.class, 
        #                             levels=spf.ordered),
         superfine.cell.class=factor(superfine.cell.class, 
                                     levels=c("ExcN","InhN","Astro","Oligo","OPC","Micro/Macro/T","Vascular","CSF"),
                                     labels=c("ExcN","InhN","Astro","Oligo","OPC","Immuno","Vasc","CSF")),
         type="snRNAseq")

seed3.94 = as.matrix(colData(spe)[,c("nmf94","nmf48")])
seed3.94 = seed3.94>0
d3.94 = cbind.data.frame(domain=as.data.frame(colData(spe))[,"domain"],
                      seed3.94) %>% 
  group_by(domain) %>% add_tally(name="total") %>%
  group_by(domain, total) %>%
  summarise_at(c("nmf94","nmf48"), sum) %>%
  tidyr::pivot_longer(c("nmf94","nmf48"), values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)

seed4.94 = as.matrix(colData(spe)[,c("nmf94","nmf48")])
seed4.94 = apply(seed4.94, 2, scale)
d4.94 = cbind.data.frame(domain=as.data.frame(colData(spe))[,"domain"],
                      seed4.94) %>% 
  group_by(domain) %>%
  summarise_at(c("nmf94","nmf48"), mean) %>% 
  tidyr::pivot_longer(c("nmf94","nmf48"), values_to="scaled.avg", names_to="nmf")

spe.94.df = left_join(d3.94[,c("domain","nmf","prop")], 
                      d4.94[,c("domain","nmf","scaled.avg")]) %>%
  mutate(type="SRT")

dot.df.94 = bind_rows(rename(sce.94.df, domain.type=superfine.cell.class), 
                      rename(spe.94.df, domain.type=domain)) %>%
  mutate(nmf_f= factor(nmf, levels=c("nmf94","nmf48")))

ggplot(dot.df.94, aes(x=nmf_f, y=domain.type, size=prop, color=scaled.avg))+
  geom_count()+theme_bw()+
  facet_wrap(vars(type), scales="free_y")+
  scale_size(range=c(0,3))+scale_color_viridis_c(option="F", direction=-1)#+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

#umaps
plotReducedDim(sce, dimred="UMAP", color_by="nmf94", point_size=.05)+
  scale_color_viridis_c(option="F", direction=-1)+
  ggtitle("nmf94")+theme(plot.title=element_text(size=18))

plotReducedDim(sce, dimred="UMAP", color_by="nmf48", point_size=.1)+
  scale_color_viridis_c(option="F", direction=-1)+
  ggtitle("nmf48")+theme(plot.title=element_text(size=18))

plotReducedDim(sce, dimred="UMAP", color_by="TTR", point_size=.1)+
  scale_color_viridis_c(option="F", direction=-1)+
  ggtitle("TTR")+theme(plot.title=element_text(size=18))


#donor bias in CP origin
brnum.palette = c("#313695","#33a02c","#fb9a99","#ff7f00","#6a3d9a",
                  "grey","#ffff99","#8dd3c7","#b2182b","grey40")
names(brnum.palette) <- c("Br2720","Br2743","Br3942","Br6423","Br6432",
                          "Br6471","Br6522","Br8325","Br8492","Br8667")

filter(as.data.frame(colData(sce)), mid.cell.class=="CSF") %>% pull(brnum) %>% table()
ggplot(filter(as.data.frame(colData(sce)), mid.cell.class=="CSF"),
       aes(x="", fill=brnum))+
  geom_bar(stat="count", position="fill", color="black")+
  scale_fill_manual(values=brnum.palette)+labs(fill="Donor", title="snRNAseq")+
  coord_polar(theta="y")+theme_void()+
  theme(plot.title=element_text(size=18), aspect.ratio=.75)

filter(as.data.frame(colData(spe)), domain=="Choroid") %>% pull(brnum) %>% table()
ggplot(filter(as.data.frame(colData(spe)), domain=="Choroid"),
       aes(x="", fill= brnum))+
  geom_bar(stat="count", position="fill", color="black")+
  scale_fill_manual(values=brnum.palette)+labs(fill="Donor", title="SRT")+
  coord_polar(theta="y")+theme_void()+
  theme(plot.title=element_text(size=18), aspect.ratio=.75)

#nmf weights by donor and CP
ggplot(mutate(as.data.frame(colData(sce)), is_CSF=mid.cell.class=="CSF"),
       aes(x=brnum, y=nmf94, fill=is_CSF))+
  geom_boxplot(outlier.size=.1)+
  theme_bw()+
  scale_fill_manual(values=c("grey50","red"))+
  labs(y="nuclei weights", x="", title="nmf94 (snRNAseq)")
ggplot(mutate(as.data.frame(colData(sce)), is_CSF=mid.cell.class=="CSF"),
       aes(x=brnum, y=nmf48, fill=is_CSF))+
  geom_boxplot(outlier.size=.1)+
  theme_bw()+
  scale_fill_manual(values=c("grey50","red"))+
  labs(y="nuclei weights", x="", title="nmf48 (snRNAseq)")

ggplot(mutate(as.data.frame(colData(spe)), is_CSF=domain=="Choroid"),
       aes(x=factor(brnum, levels=levels(spe$brnum)[c(10,1:9)]), y=nmf94, fill=is_CSF))+
  #geom_violin(scale="width")+
  geom_boxplot(outlier.size=.1)+
  theme_bw()+
  scale_fill_manual(values=c("grey50","red"))+
  labs(y="spot weights", x="", title="nmf94 (SRT)")
ggplot(mutate(as.data.frame(colData(spe)), is_CSF=domain=="Choroid"),
       aes(x=factor(brnum, levels=levels(spe$brnum)[c(10,1:9)]), y=nmf48, fill=is_CSF))+
  #geom_violin(scale="width")+
  geom_boxplot(outlier.size=.1)+
  theme_bw()+
  scale_fill_manual(values=c("grey50","red"))+
  labs(y="spot weights", x="", title="nmf48 (SRT)")


#nmf weights and TTR expr dotplot
sce$brnum_csf = ifelse(sce$mid.cell.class=="CSF", paste(sce$brnum, "CP", sep="_"),
                       paste(sce$brnum, "other", sep="_"))
sce$brnum_csf = factor(sce$brnum_csf, #levels=rev(paste0(rep(levels(sce$brnum), each=2), rep(c("_CP","_other"), length.out=20)))
                       levels=c(paste0(levels(sce$brnum), "_CP"), paste0(levels(sce$brnum), "_other"))
)

seed1 = cbind(as.matrix(colData(sce)[,c("nmf48","nmf94")]), "TTR"=logcounts(sce)["TTR",])
seed1 = seed1>0
d1 = cbind.data.frame(brnum_csf=as.data.frame(colData(sce))[,"brnum_csf"],
                      seed1) %>% 
  group_by(brnum_csf) %>% add_tally(name="total") %>%
  group_by(brnum_csf, total) %>%
  summarise_at(c("nmf48","nmf94","TTR"), sum) %>%
  tidyr::pivot_longer(c("nmf48","nmf94","TTR"), values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)

seed2 = as.matrix(colData(sce)[,c("nmf48","nmf94")])
seed2 = apply(seed2, 2, scale)
seed2 = cbind(seed2, "TTR"=logcounts(sce)["TTR",])
d2 = cbind.data.frame(brnum_csf=as.data.frame(colData(sce))[,"brnum_csf"],
                      seed2) %>% 
  group_by(brnum_csf) %>%
  summarise_at(c("nmf48","nmf94","TTR"), mean) %>% 
  tidyr::pivot_longer(c("nmf48","nmf94","TTR"), values_to="scaled.avg", names_to="nmf")

dot.df = left_join(d1[,c("brnum_csf","nmf","n")], 
                   d2[,c("brnum_csf","nmf","scaled.avg")]) %>%
  mutate(brnum_csf=factor(brnum_csf, levels=levels(sce$brnum_csf)),
         y_labs=factor(nmf, levels=c("TTR","nmf94","nmf48")))

nmf.plot <- ggplot(dot.df, aes(x=brnum_csf, y=y_labs, color=scaled.avg, size=n))+
  geom_count()+theme_bw()+
  scale_size(range=c(0,6))+#, limits=c(0,1), breaks=c(0,.5,1))+
  labs(title="snRNA-seq")+
  scale_color_viridis_c(option="F", direction=-1, limits=c(min(filter(dot.df, nmf!="TTR")$scaled.avg),
                                                           max(filter(dot.df, nmf!="TTR")$scaled.avg)))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), legend.key.size=unit(12,"pt"))

expr.plot <- ggplot(dot.df, aes(x=brnum_csf, y=y_labs, color=scaled.avg, size=n))+
  geom_count()+theme_bw()+
  scale_size(range=c(0,6))+#, limits=c(0,1), breaks=c(0,.5,1))+
  labs(title="snRNA-seq")+
  scale_color_gradient(low="grey90", high="black", limits=c(min(filter(dot.df, nmf=="TTR")$scaled.avg),
                                                           max(filter(dot.df, nmf=="TTR")$scaled.avg)))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), legend.key.size=unit(12,"pt"),
        axis.text.y=element_text(face="italic"))

pdf(file="plots/revision/supp_nmf94-donor-bias_dotplot-nmf-TTR.pdf", height=6, width=8)
gridExtra::grid.arrange(nmf.plot, expr.plot, ncol=1)
dev.off()




spe$brnum_csf = ifelse(spe$domain=="Choroid", paste(spe$brnum, "CP", sep="_"),
                       paste(spe$brnum, "other", sep="_"))
table(spe$brnum_csf)
spe$brnum_csf = factor(spe$brnum_csf, #levels=rev(paste0(rep(levels(sce$brnum), each=2), rep(c("_CP","_other"), length.out=20)))
                       levels=levels(sce$brnum_csf))

seed3 = cbind(as.matrix(colData(spe)[,c("nmf48","nmf94")]), "TTR"=logcounts(spe)["TTR",])
seed3 = seed3>0
d3 = cbind.data.frame(brnum_csf=as.data.frame(colData(spe))[,"brnum_csf"],
                      seed3) %>% 
  group_by(brnum_csf) %>% add_tally(name="total") %>%
  group_by(brnum_csf, total) %>%
  summarise_at(c("nmf48","nmf94","TTR"), sum) %>%
  tidyr::pivot_longer(c("nmf48","nmf94","TTR"), values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)

d3 = bind_rows(d3, filter(d3, brnum_csf=="Br2720_CP") %>% mutate(n=0, prop=0, brnum_csf="Br6432_CP"))

seed4 = as.matrix(colData(spe)[,c("nmf48","nmf94")])
seed4 = apply(seed4, 2, scale)
seed4 = cbind(seed4, "TTR"=logcounts(spe)["TTR",])
d4 = cbind.data.frame(brnum_csf=as.data.frame(colData(spe))[,"brnum_csf"],
                      seed4) %>% 
  group_by(brnum_csf) %>%
  summarise_at(c("nmf48","nmf94","TTR"), mean) %>% 
  tidyr::pivot_longer(c("nmf48","nmf94","TTR"), values_to="scaled.avg", names_to="nmf")

d4 = bind_rows(d4, filter(d4, brnum_csf=="Br2720_CP") %>% mutate(scaled.avg=0, brnum_csf="Br6432_CP"))


dot.df2 = left_join(d3[,c("brnum_csf","nmf","n","prop")], 
                   d4[,c("brnum_csf","nmf","scaled.avg")]) %>%
  mutate(brnum_csf=factor(brnum_csf, levels=levels(spe$brnum_csf)),
         y_labs=factor(nmf, levels=c("TTR","nmf94","nmf48")))

nmf.plot2 <- ggplot(dot.df2, aes(x=brnum_csf, y=y_labs, color=scaled.avg, size=prop))+
  geom_count()+theme_bw()+#geom_text(aes(label=n), color="white")+
  scale_size(range=c(0,6), limits=c(0,1), breaks=c(0,.5,1))+
  labs(title="SRT")+
  scale_color_viridis_c(option="F", direction=-1, limits=c(min(filter(dot.df2, nmf!="TTR")$scaled.avg),
                                                           max(filter(dot.df2, nmf!="TTR")$scaled.avg)))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), legend.key.size=unit(12,"pt"))

expr.plot2 <- ggplot(dot.df2, aes(x=brnum_csf, y=y_labs, color=scaled.avg, size=prop))+
  geom_count()+theme_bw()+#geom_text(aes(label=n), color="white")+
  scale_size(range=c(0,6), limits=c(0,1), breaks=c(0,.5,1))+
  labs(title="SRT")+
  scale_color_gradient(low="grey90", high="black", limits=c(min(filter(dot.df2, nmf=="TTR")$scaled.avg),
                                                            max(filter(dot.df2, nmf=="TTR")$scaled.avg)))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), legend.key.size=unit(12,"pt"),
        axis.text.y=element_text(face="italic"))

pdf(file="plots/revision/supp_nmf94-donor-bias_SRT-dotplot-nmf-TTR.pdf", height=6, width=8)
gridExtra::grid.arrange(nmf.plot2, expr.plot2, ncol=1)
dev.off()

#spot plots
sub=spe[,spe$brnum=="Br3942"]
sub$is_choroid = sub$domain=="Choroid"

plotVisium(sub, annotate = "TTR", highlight = "is_choroid", facets = "sample_id",
           assay = "logcounts", image = FALSE, #x_coord = "array_row", y_coord= "array_col"
           )+
  scale_color_manual(values=c("grey","#00006a"))+
  scale_fill_gradient(low="grey90", high="black")+ggtitle("TTR")+
  theme(plot.title=element_text(size=18, hjust=.5), strip.text=element_blank())
plotVisium(sub, annotate = "nmf94", highlight = "is_choroid", facets = "sample_id", image = FALSE)+
  scale_color_manual(values=c("grey","#00006a"))+
  scale_fill_gradient(low="grey90", high="black")+ggtitle("nmf94")+
  theme(plot.title=element_text(size=18, hjust=.5), strip.text=element_blank())
plotVisium(sub, annotate = "nmf48", highlight = "is_choroid", facets = "sample_id", image = FALSE)+
  scale_color_manual(values=c("grey","#00006a"))+
  scale_fill_gradient(low="grey90", high="black")+ggtitle("nmf48")+
  theme(plot.title=element_text(size=18, hjust=.5), strip.text=element_blank())
