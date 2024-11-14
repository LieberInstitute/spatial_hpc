library(scater)
library(BiocSingular)
library(scry)
library(dplyr)
library(ggplot2)
library(here)

#i couldn't reproduce the code found here:
#https://github.com/LieberInstitute/spatial_hpc/blob/926d84cca64335bd1925600c7855a61c0de3b041/snRNAseq_hpc/code/build_sce/feature_selection_dimred_clustering.R

#i couldn't reproduce his work but thank god he actually saved the object before subsetting it!!!!
load("snRNAseq_hpc/processed-data/sce/sce_clustered_round2.rda")

df = as.data.frame(colData(sce)[,c(1:20,23,24)]) %>% mutate(syt1= logcounts(sce)["SYT1",])
ggplot(mutate(df, nrn_clus=paste(neuron, k_5_louvain_initial)), 
       aes(x=nrn_clus, y=syt1))+
  geom_boxplot(outlier.size=.5)+theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1))

ordered_clus = filter(df, neuron==F) %>% group_by(k_5_louvain_initial) %>% summarise(med_det=median(detected)) %>% arrange(desc(med_det))
ggplot(filter(df, neuron==F) %>%
         mutate(k_5_louvain_initial=factor(k_5_louvain_initial, levels=ordered_clus$k_5_louvain_initial),
                low_nrn = if_else(k_5_louvain_initial %in% c("6","8","46"), as.character(k_5_louvain_initial), "none")), 
       aes(x=k_5_louvain_initial, y=detected, fill=low_nrn))+
  geom_boxplot(outlier.size=.5)+theme_bw()+
  scale_fill_manual(values=c("red3","darkblue","violet","grey"))

sce$low_nrn = if_else(sce$k_5_louvain_initial %in% c("6","8","46"), as.character(sce$k_5_louvain_initial), "none")
p1 <- plotReducedDim(sce, dimred="UMAP", color_by="low_nrn", point_size=.3)+scale_color_manual(values=c("red3","darkblue","violet","grey"))
p2 <- plotReducedDim(sce, dimred="UMAP", color_by="detected", point_size=.3)
gridExtra::grid.arrange(p1, p2, ncol=2)

tmp = sce[,sce$neuron==TRUE]
tmp$k_5_louvain_initial = droplevels(tmp$k_5_louvain_initial)
tmp$k_5_louvain_initial = factor(tmp$k_5_louvain_initial, levels=c(38,59,
                                                                   20,40,34,
                                                                   43,
                                                                   9,21,31,
                                                                   42,1,10,26,
                                                                   28,
                                                                   2,
                                                                   37
                                                                   ))

plotDots(tmp, features=c("PRLR","ABCB1","FLT1","ATP10A","AKAP12","CFAP73","GLIS3","RFX4","MBP","CTNNA3","LHFPL3","TBXAS1"), 
         group="k_5_louvain_initial")+
  scale_color_viridis_c(option="F", direction=-1)

#final round
load("snRNAseq_hpc/processed-data/sce/sce_clustered_round4.rda")
df2 = as.data.frame(colData(sce)[,c(1:20,23,24,28)]) %>% mutate(syt1= logcounts(sce)["SYT1",])
ggplot(mutate(df2, nrn_clus=paste(sapply(neuron, isFALSE), k_5_louvain)), 
       aes(x=nrn_clus, y=syt1))+
  geom_boxplot(outlier.size=.5)+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1))

ordered_clus2 = filter(df2, neuron==T) %>% group_by(k_5_louvain) %>% summarise(med_det=median(detected)) %>% arrange(desc(med_det))
ggplot(filter(df2, neuron==T) %>% 
         mutate(k_5_louvain=factor(k_5_louvain, levels=ordered_clus2$k_5_louvain),
                low_nrn = if_else(k_5_louvain %in% c("17"), as.character(k_5_louvain), "none")),# %>%
         #mutate(low_nrn = if_else(k_5_louvain_initial %in% c("6","8","46"), as.character(k_5_louvain_initial), "none")), 
       aes(x=k_5_louvain, y=detected, fill=low_nrn))+
  geom_boxplot(outlier.size=.5)+theme_bw()+
  scale_fill_manual(values=c("red3","grey"))

tmp = sce[,sce$neuron==FALSE]
tmp$k_5_louvain = droplevels(tmp$k_5_louvain)
tmp$k_5_louvain = factor(tmp$k_5_louvain, levels=c(34,36,62,
                                                    19,47,41,30,
                                                    48,
                                                    8,20,27,
                                                    1,9,
                                                    44,24,
                                                    2,6,42,
                                                    23
))
plotDots(tmp, features=c("PRLR","ABCB1","FLT1","ATP10A","AKAP12","CFAP73","GLIS3","RFX4","MBP","CTNNA3","LHFPL3","TBXAS1"), 
         group="k_5_louvain")+
  scale_color_viridis_c(option="F", direction=-1)

setdiff(unique(tmp$k_5_louvain), c(34,36,62,
  19,47,41,30,
  48,
  8,20,27,
  1,9,
  44,24,
  2,6,42,
  23
))
