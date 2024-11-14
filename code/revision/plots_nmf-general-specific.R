library(SingleCellExperiment)
library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(scater)

set.seed(123)

#sce
load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))
#baseline non-zero
non0.nuc = colSums(as.matrix(colData(sce)[,paste0("nmf",1:100)])>0)
plot(ecdf(log10(non0.nuc)), xlim=c(0,5), xlab="# non-zero weighted nuclei per NMF")

#spe
load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))
non0.spots = colSums(as.matrix(colData(spe)[,paste0("nmf",1:100)])>0)
plot.spots = ifelse(na.exclude(non0.spots<1050), "red", "black")
plot(ecdf(log10(non0.spots)), xlim=c(0,5), xlab="# non-zero weighted spots per NMF",
     col=sort(plot.spots, decreasing=T), 
     main="NMF with <1050 spots removed")

remove.nmf = names(non0.spots[non0.spots<1050])
remove.nmf[c(2,3,6)] <- c("nmf2","nmf3","nmf16")
remove.nmf = c(remove.nmf, "nmf37", "nmf28")

##############
## dotplot

seed1 = as.matrix(colData(sce)[,paste0("nmf",1:100)])
seed1 = seed1>0
d1 = cbind.data.frame(superfine.cell.class=as.data.frame(colData(sce))[,"superfine.cell.class"],
                      seed1) %>% 
  group_by(superfine.cell.class) %>% add_tally(name="total") %>%
    group_by(superfine.cell.class, total) %>%
  summarise_at(paste0("nmf",1:100), sum) %>%
  tidyr::pivot_longer(paste0("nmf",1:100), values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)

seed2 = as.matrix(colData(sce)[,paste0("nmf",1:100)])
seed2 = apply(seed2, 2, scale)
d2 = cbind.data.frame(superfine.cell.class=as.data.frame(colData(sce))[,"superfine.cell.class"],
                      seed2) %>% 
  group_by(superfine.cell.class) %>%
  summarise_at(paste0("nmf",1:100), mean) %>% 
  tidyr::pivot_longer(paste0("nmf",1:100), values_to="scaled.avg", names_to="nmf")


spf.ordered=c("GC.3","GC.1","GC.2","GC.4","GC.5",
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
         "Ependy","CP.1","CP.2","CP.3","Endo.2","Endo.1","PC/SMC","VLMC")

dot.df = left_join(d1[,c("superfine.cell.class","nmf","prop")], 
                   d2[,c("superfine.cell.class","nmf","scaled.avg")]) %>%
  mutate(superfine.cell.class=factor(superfine.cell.class, 
                                     levels=spf.ordered))

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


dot.df$nmf_f = factor(dot.df$nmf, levels=c(nmf.ordered.remove, nmf.ordered.keep))

ggplot(dot.df, aes(x=nmf_f, y=superfine.cell.class, size=prop, color=scaled.avg))+
  geom_count()+theme_bw()+
  scale_size(range=c(0,3))+scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

#possible addition of correlation with library size, etc.
libsize.cor = sapply(nmf.ordered, function(x) cor(sce$sum,colData(sce)[,x]))
names(libsize.cor) <- nmf.ordered
detected.cor = sapply(nmf.ordered, function(x) cor(sce$detected,colData(sce)[,x]))
names(detected.cor) <- nmf.ordered
mito.cor = sapply(nmf.ordered, function(x) cor(sce$subsets_Mito_percent,colData(sce)[,x]))
names(mito.cor) <- nmf.ordered

dot.df$new.y = as.character(dot.df$superfine.cell.class)
add.y = filter(dot.df, superfine.cell.class=="GC.1") %>% mutate(new.y="GC.1")
add.y1 = left_join(filter(dot.df, superfine.cell.class=="GC.1") %>% mutate(new.y="libsize.cor", prop=1) %>%
                     select(superfine.cell.class, nmf, prop, nmf_f, new.y),
                   cbind.data.frame("nmf"=names(libsize.cor), "scaled.avg"=libsize.cor))
add.y2 = left_join(filter(dot.df, superfine.cell.class=="GC.1") %>% mutate(new.y="detected.cor", prop=1) %>%
                     select(superfine.cell.class, nmf, prop, nmf_f, new.y),
                   cbind.data.frame("nmf"=names(detected.cor), "scaled.avg"=detected.cor))
add.y3 = left_join(filter(dot.df, superfine.cell.class=="GC.1") %>% mutate(new.y="mito.cor", prop=1) %>%
                     select(superfine.cell.class, nmf, prop, nmf_f, new.y),
                   cbind.data.frame("nmf"=names(mito.cor), "scaled.avg"=mito.cor))

ggplot(bind_rows(add.y, add.y1, add.y2, add.y3), 
       aes(x=nmf_f, y=new.y, size=prop, color=scaled.avg))+
  geom_count()+theme_bw()+
  scale_size(range=c(0,3))+
  #scale_color_gradientn(colors=RColorBrewer::brewer.pal(n=7, "RdYlBu")[7:1], limits=c(-.8,.8))+
  scale_color_gradient2(low="#4575b4", mid="grey90", high="#d73027", limits=c(-.4,.8))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))


#spe dotplot
seed3 = as.matrix(colData(spe)[,nmf.ordered.keep])
seed3 = seed3>0
d3 = cbind.data.frame(domain=as.data.frame(colData(spe))[,"domain"],
                      seed3) %>% 
  group_by(domain) %>% add_tally(name="total") %>%
  group_by(domain, total) %>%
  summarise_at(nmf.ordered.keep, sum) %>%
  tidyr::pivot_longer(nmf.ordered.keep, values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)

seed4 = as.matrix(colData(spe)[,nmf.ordered.keep])
seed4 = apply(seed4, 2, scale)
d4 = cbind.data.frame(domain=as.data.frame(colData(spe))[,"domain"],
                      seed4) %>% 
  group_by(domain) %>%
  summarise_at(nmf.ordered.keep, mean) %>% 
  tidyr::pivot_longer(nmf.ordered.keep, values_to="scaled.avg", names_to="nmf")


dot.df2 = left_join(d3[,c("domain","nmf","prop")], 
                    d4[,c("domain","nmf","scaled.avg")])

dot.df2$nmf_f = factor(dot.df2$nmf, levels=nmf.ordered.keep)

ggplot(dot.df2, aes(x=nmf_f, y=domain, size=prop, color=scaled.avg))+
  geom_count()+theme_bw()+
  scale_size(range=c(0,3))+scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

#### add corr of lib size
libsize.cor2 = sapply(nmf.ordered.keep, function(x) cor(spe$sum_umi,colData(spe)[,x]))
names(libsize.cor2) <- nmf.ordered.keep
detected.cor2 = sapply(nmf.ordered.keep, function(x) cor(spe$sum_gene,colData(spe)[,x]))
names(detected.cor2) <- nmf.ordered.keep
mito.cor2 = sapply(nmf.ordered.keep, function(x) cor(spe$expr_chrM_ratio,colData(spe)[,x]))
names(mito.cor2) <- nmf.ordered.keep

dot.df2$new.y = as.character(dot.df2$domain)
add.y = filter(dot.df2, domain=="GCL") %>% mutate(new.y="GCL")
add.y1 = left_join(filter(dot.df2, domain=="GCL") %>% mutate(new.y="libsize.cor", prop=1) %>%
                     select(domain, nmf, prop, nmf_f, new.y),
                   cbind.data.frame("nmf"=names(libsize.cor2), "scaled.avg"=libsize.cor2))
add.y2 = left_join(filter(dot.df2, domain=="GCL") %>% mutate(new.y="detected.cor", prop=1) %>%
                     select(domain, nmf, prop, nmf_f, new.y),
                   cbind.data.frame("nmf"=names(detected.cor2), "scaled.avg"=detected.cor2))
add.y3 = left_join(filter(dot.df2, domain=="GCL") %>% mutate(new.y="mito.cor", prop=1) %>%
                     select(domain, nmf, prop, nmf_f, new.y),
                   cbind.data.frame("nmf"=names(mito.cor2), "scaled.avg"=mito.cor2))

ggplot(bind_rows(add.y, add.y1, add.y2, add.y3), 
       aes(x=nmf_f, y=new.y, size=prop, color=scaled.avg))+
  geom_count()+theme_bw()+
  scale_size(range=c(0,3))+
  #scale_color_gradientn(colors=RColorBrewer::brewer.pal(n=7, "RdYlBu")[7:1], limits=c(-.8,.8))+
  scale_color_gradient2(low="#4575b4", mid="grey90", high="#d73027", limits=c(-.6,1))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
