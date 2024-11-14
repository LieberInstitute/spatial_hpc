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
test$xlabel = paste(test$brnum, test$round, test$sort)

ggplot(group_by(test, nmf_f, xlabel) %>% summarise(n.non0.spots=sum(n.non0.spots)), 
       aes(y=n.non0.spots, x=nmf_f, fill=xlabel))+
  geom_bar(stat="identity", position="fill", color="grey50")+theme_bw()+
  scale_fill_manual(values=c("#a6cee3","#74add1","#4575b4","#313695",
                             "#b2df8a","#33a02c",
                             "#fb9a99","#e31a1c",#3492
                             "#fee08b","#fdbf6f","#ff7f00","#f46d43",#6423
                             "#c2a5cf","#6a3d9a",#6432
                             "grey","grey50",#6471
                             "#ffff99","goldenrod","#b35806","#7f3b08",#6522
                             "#8dd3c7","aquamarine",#8325
                             "#f4a582","#b2182b",#8492
                             "grey40","black"))+#8667
  labs(y="% nuclei with non-0 loading", title="General NMF patterns",fill="Sample")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

#specific patterns only
spec.nmf = setdiff(nmf.ordered.keep, gen.nmf)
spec.non0.mtx = as.matrix(colData(sce)[,spec.nmf])
spec.non0.mtx = spec.non0.mtx>0
spec.non0 = cbind(as.data.frame(spec.non0.mtx), 
                  as.data.frame(colData(sce)[rownames(spec.non0.mtx),c("brnum","round","sort")]))
test2 = group_by(spec.non0, brnum, round, sort) %>% summarise_at(spec.nmf, sum) %>%
  tidyr::pivot_longer(all_of(spec.nmf), names_to="nmf", values_to="n.non0.spots")
test2$nmf_f = factor(test2$nmf, levels=nmf.ordered.keep)
test2$xlabel = paste(test2$brnum, test2$round, test2$sort)

ggplot(group_by(test2, nmf_f, xlabel) %>% summarise(n.non0.spots=sum(n.non0.spots)), 
       aes(y=n.non0.spots, x=nmf_f, fill=xlabel))+
  geom_bar(stat="identity", position="fill", color="grey50")+theme_bw()+
  scale_fill_manual(values=c("#a6cee3","#74add1","#4575b4","#313695",#2720
                             "#b2df8a","#33a02c",#2743
                             "#fb9a99","#e31a1c",#3492
                             "#fdbf6f","#ff7f00","#fee08b","#f46d43",#6423
                             "#c2a5cf","#6a3d9a",#6432
                             "grey","grey50",#6471
                             "#ffff99","goldenrod","#b35806","#7f3b08",#6522
                             "#8dd3c7","aquamarine",#8325
                             "#f4a582","#b2182b",#8492
                             "grey40","black"))+#8667
  labs(y="% nuclei with non-0 loading", title="Specific NMF patterns",fill="Sample")

#################
## ok now do same as above but with SRT
#################

# general patterns
gen.non0.mtx = as.matrix(colData(spe)[,gen.nmf])
gen.non0.mtx = gen.non0.mtx>0
gen.non0 = cbind(as.data.frame(gen.non0.mtx), 
                 as.data.frame(colData(spe)[rownames(gen.non0.mtx),c("brnum","sample_id")]))
test = group_by(gen.non0, brnum, sample_id) %>% summarise_at(gen.nmf, sum) %>%
  tidyr::pivot_longer(all_of(gen.nmf), names_to="nmf", values_to="n.non0.spots")
test$nmf_f = factor(test$nmf, levels=nmf.ordered.keep)
test$brnum_f = factor(test$brnum, levels=levels(spe$brnum)[c(10,1:9)])

ggplot(group_by(test, nmf_f, brnum_f) %>% summarise(n.non0.spots=sum(n.non0.spots)), 
       aes(x=n.non0.spots, y=nmf_f, fill=brnum_f))+
  geom_bar(stat="identity", position="fill", color="grey50")+theme_bw()+
  scale_fill_manual(values=c("#313695",#2720
                             "#33a02c",#2743
                             "#fb9a99",#3492
                             "#ff7f00",#6423
                             "#6a3d9a",#6432
                             "grey",#6471
                             "#ffff99",#6522
                             "#8dd3c7",#8325
                             "#b2182b",#8492
                             "grey40","black"#8667
                             ))+
  labs(x="% spots with non-0 loading", title="General NMF patterns",fill="Donor")

#specific patterns
spec.nmf = setdiff(nmf.ordered.keep, gen.nmf)
spec.non0.mtx = as.matrix(colData(spe)[,spec.nmf])
spec.non0.mtx = spec.non0.mtx>0
spec.non0 = cbind(as.data.frame(spec.non0.mtx), 
                  as.data.frame(colData(spe)[rownames(spec.non0.mtx),c("brnum","sample_id")]))
test2 = group_by(spec.non0, brnum) %>% summarise_at(spec.nmf, sum) %>%
  tidyr::pivot_longer(all_of(spec.nmf), names_to="nmf", values_to="n.non0.spots")
test2$nmf_f = factor(test2$nmf, levels=nmf.ordered.keep)
test2$brnum_f = factor(test2$brnum, levels=levels(spe$brnum)[c(10,1:9)])

ggplot(group_by(test2, nmf_f, brnum_f) %>% summarise(n.non0.spots=sum(n.non0.spots)), 
       aes(x=n.non0.spots, y=nmf_f, fill=brnum_f))+
  geom_bar(stat="identity", position="fill", color="grey50")+theme_bw()+
  scale_fill_manual(values=c("#313695",#2720
                             "#33a02c",#2743
                             "#fb9a99",#3492
                             "#ff7f00",#6423
                             "#6a3d9a",#6432
                             "grey",#6471
                             "#ffff99",#6522
                             "#8dd3c7",#8325
                             "#b2182b",#8492
                             "grey40","black"#8667
  ))+
  labs(x="% spots with non-0 loading", title="Specific NMF patterns",fill="Donor")

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
ggplot(filter(as.data.frame(colData(sce)), mid.cell.class=="CSF"),
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
  ))+labs(fill="Donor", title="snRNAseq")+
  coord_polar(theta="y")+theme_void()+
  theme(plot.title=element_text(size=18), aspect.ratio=1)

ggplot(filter(as.data.frame(colData(spe)), domain=="Choroid"),
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
  ))+labs(fill="Donor", title="SRT")+
  coord_polar(theta="y")+theme_void()+
  theme(plot.title=element_text(size=18), aspect.ratio=1)

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
