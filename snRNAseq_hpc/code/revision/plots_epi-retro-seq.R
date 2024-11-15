library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scater)
set.seed(123)

load("snRNAseq_hpc/processed-data/revision/sce_epiretroseq_filtered.rda")

table(mch$Source)
table(mch$Cocluster)
table(mch$Subclass)
mch$Subclass = droplevels(mch$Subclass)
table(colData(mch)[,c("Source","Subclass")])
colnames(colData(mch))


colSums(as.matrix(colData(mch)[,paste0("nmf",1:100)]))
colSums(!is.na(colData(mch2)[,paste0("nmf",1:100)]))

mch = mch[,mch$Source!="DGa" & mch$Source!="DGp"]

nmf.use = c("nmf52","nmf11","nmf63","nmf61",#mc, ca3, ca2
  "nmf15",#ca1
  "nmf32","nmf40","nmf54",#pros/sub
  "nmf65",#deep sub
  "nmf53","nmf22", #L6s
  "nmf51","nmf68",#l5
  "nmf27","nmf45","nmf84","nmf17"#L2/3
)

non0.nuc = colSums(as.matrix(colData(mch)[,paste0("nmf",1:100)])>0)

plot1 = ifelse(na.exclude(names(non0.nuc) %in% nmf.use), "red", "black")
plot(ecdf(non0.spots), xlab="# non-zero weighted spots per NMF", cex=1,
     col=plot1,#sort(plot.spots, decreasing=T), 
     main="Red: CA, sub, RHP patterns")

non0.use = non0.spots[nmf.use]
plot.spots = ifelse(na.exclude(non0.use<45), "red", "black")
plot(ecdf(non0.use), xlab="# non-zero weighted spots per NMF",
     col=sort(plot.spots, decreasing=T), 
     main="NMF with <45 spots removed")

non0.keep = names(non0.use)[non0.use>45]

#tissue source dotplot
seed1 = as.matrix(colData(mch)[,non0.keep])
seed1 = seed1>0
d1 = cbind.data.frame(as.data.frame(colData(mch))[,c("Source","Target")],
                      seed1) %>% 
  group_by(Source) %>% add_tally(name="total") %>%
  group_by(Source, total) %>%
  summarise_at(non0.keep, sum) %>%
  tidyr::pivot_longer(non0.keep, values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)

seed2 = as.matrix(colData(mch)[,non0.keep])
seed2 = apply(seed2, 2, scale)
d2 = cbind.data.frame(as.data.frame(colData(mch))[,c("Source","Target")],
                      seed2) %>% 
  group_by(Source) %>%
  summarise_at(non0.keep, mean) %>% 
  tidyr::pivot_longer(non0.keep, values_to="scaled.avg", names_to="nmf")

nmf.ordered = factor(non0.keep, levels=non0.keep)

dot.df = left_join(d1[,c("Source","nmf","prop","n")], 
                   d2[,c("Source","nmf","scaled.avg")]) %>%
  mutate(nmf_f=factor(nmf, levels=nmf.ordered),
         Source_f=factor(Source, levels=c("DGp","DGa","CAp","CAa","ENT")))

ggplot(dot.df, aes(x=nmf_f, y=Source_f, size=n, color=scaled.avg))+
  geom_count()+theme_bw()+
  labs(y="Collection source", x="CA & periallocortex NMF patterns", size="n nuclei",
       title="mouse epi-retroSeq")+
  scale_size(range=c(0,6))+scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        text=element_text(size=16))


#target region dotplot
d3 = cbind.data.frame(as.data.frame(colData(mch))[,c("Source","Target")],
                      seed1) %>% 
  filter(Source!="DGa", Source!="DGp") %>%
  group_by(Target) %>% add_tally(name="total") %>%
  group_by(Target, total) %>%
  summarise_at(non0.keep, sum) %>%
  tidyr::pivot_longer(non0.keep, values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)
d4 = cbind.data.frame(as.data.frame(colData(mch))[,c("Source","Target")],
                      seed2) %>% 
  filter(Source!="DGa", Source!="DGp") %>%
  group_by(Target) %>%
  summarise_at(non0.keep, mean) %>% 
  tidyr::pivot_longer(non0.keep, values_to="scaled.avg", names_to="nmf")

dot.df2 = left_join(d3[,c("Target","nmf","prop","n")], 
                   d4[,c("Target","nmf","scaled.avg")]) %>%
  mutate(nmf_f=factor(nmf, levels=nmf.ordered))#,
         #Source_f=factor(Source, levels=c("DGp","DGa","CAp","CAa","ENT")))

dot.df2$Target_f = factor(dot.df2$Target, levels=c("HPF","ENT","MOB",
                                                   "RSP","PTLp","ACA",
                                                   "PIR","STR","TH","AMY","HY",
                                                   "MOp","PFC"
                                                   ))

ggplot(dot.df2, aes(x=nmf_f, y=Target_f, size=n, color=scaled.avg))+
  geom_count()+theme_bw()+labs(y="Target region", x="CA & periallocortex NMF patterns", size="# nuclei",
                               title="mouse epi-retroSeq")+
  scale_size(range=c(0,6))+scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
        text=element_text(size=16))
