library(SpatialExperiment)
library(ggspavis)
library(dplyr)
library(ggplot2)


################## BayesSpace
#can start with  bayesspace spe because it has the precast cluster results in it
load("processed-data/06_clustering/BayesSpace/revision/spe_bayes_k18-kmeans-10k.Rdata")
spe$amy<-ifelse(spe$sample_id %in%
                  c('Br6423_V10B01-085_C1','Br6432_V10B01-086_B1') &
                  spe$domain %in% c('RHP'),T,F)
spe$amy<-ifelse(spe$sample_id =='Br6423_V10B01-085_A1',T,spe$amy)
spe$thal<-ifelse(spe$sample_id=='Br8325_V11A20-297_B1' & spe$domain %in% c('SUB','SUB.RHP'),T,F)

spe$cluster = as.character(spe$cluster)
spe$cluster = ifelse(spe$amy==TRUE, "Amygdala",spe$cluster)
spe$cluster = ifelse(spe$thal==TRUE, "Thalamus",spe$cluster)

spe$cluster_f = factor(spe$cluster, levels=c("Thalamus","GABA","Amygdala",
                                             "RHP","SUB.RHP","SUB",
                                             "CA1.1","CA1.2","CA2.4.1","CA2.4.2",
                                             "GCL",
                                             "ML","SR.SLM","SL.SR","SLM.SGZ",
                                             "WM.1","WM.2","WM.3","Vascular","Choroid"))

spe$spatial.cluster.f = factor(spe$spatial.cluster, levels=c(16,5,18,3,8,15,4,13,
                               14,9,11,10,1,
                               7,12,17,6,2
))


plot.df = as.data.frame(colData(spe)[,c("cluster_f", "spatial.cluster.f")]) %>% 
  group_by(domain) %>% add_tally(name="total.precast") %>%
  group_by(domain, spatial.cluster.f, total.precast, .drop=FALSE) %>% tally() %>%
  mutate(prop.precast=n/total.precast, prop.precast=ifelse(is.na(prop.precast), 0, prop.precast))

ggplot(plot.df, aes(x=spatial.cluster.f, y=domain, fill=prop.precast))+
  geom_tile()+
  scale_fill_viridis_c(option="F", limits=c(0,1), breaks=c(0,.25,.5,.75,1))+
  labs(fill="prop. PRECAST\ndomain nuclei", y="PRECAST domain", x="BayesSpace cluster")+
  theme_minimal()+theme(text=element_text(size=16), legend.title=element_text(size=12),
                        legend.text=element_text(size=10))

plot.df2 = as.data.frame(colData(spe)[,c("cluster_f", "spatial.cluster.f")]) %>% 
  group_by(cluster_f, spatial.cluster.f, .drop=FALSE) %>% 
  tally() %>% mutate(log10.n=log10(n+1))

ggplot(plot.df2, aes(x=spatial.cluster.f, y=cluster_f, fill=log10.n))+
  geom_tile()+
  scale_fill_gradient2(low="white", mid="skyblue",high="black", midpoint=2)+#median(plot.df2$log10.n))+
  #scale_fill_viridis_c(option="F")#, limits=c(0,1), breaks=c(0,.25,.5,.75,1))+
  labs(fill="log10(# spots)", y="PRECAST cluster", x="BayesSpace cluster")+
  theme_minimal()+theme(text=element_text(size=16), legend.title=element_text(size=12),
                        legend.text=element_text(size=10))


################## GraphST

graphst = read.csv("processed-data/06_clustering/GraphST/hpc_allSamples_graphst_k16.csv") %>%
  mutate(spot_id = gsub("1-[0-9]*", "1_", spot_id))
graphst$spot_id = paste0(graphst$spot_id, graphst$captureArea)
colSums(is.na(graphst))

#fix spot codes in colData
spe$spot_id = paste(sapply(strsplit(rownames(colData(spe)), "_"), function(x) x[[1]]), spe$slide, spe$array, sep="_")

matching = left_join(as.data.frame(colData(spe))[,c("spot_id","brnum","PRECAST_k18","cluster_f")] %>%
                      mutate(brnum=as.character(brnum)),
          graphst, by=c("spot_id", "brnum"="donor_id"))
colSums(is.na(matching))
colnames(matching)

l1 = colnames(full.df)[7:13]
names(l1) = l1
lapply(l1, function(x) {
  tmp = full.df[,c("domain",x)]
  colnames(tmp) = c("domain","graphst")
  group_by(tmp, graphst) %>% add_tally(name="total.cluster") %>%
    group_by(domain, graphst, total.cluster) %>% tally() %>% ungroup() %>%
    mutate(prop=n/total.cluster) %>% 
    group_by(domain) %>% slice_max(n=3, prop) %>% filter(domain=="RHP")
})

#lamb_10 has best matching with rhp but cluster_lamb_0_1_and_1 is what is in the code as used:
#https://github.com/LieberInstitute/spatial_hpc/blob/926d84cca64335bd1925600c7855a61c0de3b041/plots/06_clustering/GraphST/add_graphst_to_spe.R
summ.df = group_by(full.df, cluster_lamb_0_1_and_1) %>% add_tally(name="total.cluster") %>%
  group_by(domain, cluster_lamb_0_1_and_1, total.cluster) %>% tally() %>% ungroup() %>%
  mutate(prop=n/total.cluster) %>% 
  group_by(domain) %>% slice_max(n=3, prop)
print(n=48, summ.df)

identical(spe$spot_id, full.df$spot_id)
spe$cluster_lamb_10 = full.df$cluster_lamb_10
#lamb_10
spe$graphst.f = factor(spe$cluster_lamb_10, levels= c(14,
         2,
         12,
         3,
         13,
         16,#SUB.RHP 49%
         11,#12,#RHP 54%
         4,#GABA 43%
         #NA,#SLSR 6 = 13%
         6,7,
         #6,#SR.SLM 53%
         #NA,#SLM.SGZ
         10,8,15,
         9,#vascular
         1,#choroid
         5
         ))
#lamb_0_1_and_1
spe$cluster_lamb_0_1_and_1 = full.df$cluster_lamb_0_1_and_1
spe$graphst.f = factor(spe$cluster_lamb_0_1_and_1, levels= c(11,
                                                             9,2,
                                                             7,
                                                             3,
                                                             15,
                                                             14,#SUB.RHP
                                                             #7,#RHP
                                                             4,#GABA
                                                             8,#SL.SR
                                                             6,#SR.SLM
                                                             5,#SLM.SGZ
                                                             13,12,16,
                                                             1,10
                                                             
))

matching$graphst.f = factor(matching$cluster_lamb_0_1_and_1, levels= c(11,
                                                                 9,2,
                                                                 7,
                                                                 3,
                                                                 15,
                                                                 14,#SUB.RHP
                                                                 #7,#RHP
                                                                 4,#GABA
                                                                 8,#SL.SR
                                                                 6,#SR.SLM
                                                                 5,#SLM.SGZ
                                                                 13,12,16,
                                                                 1,10
))

matching$graphst.f = factor(matching$cluster_lamb_0_1_and_1, levels=c(4,7,14,15,3,9,2,11,
                                                                      8,6,13,5,12,
                                                                      16,10,1
))

plot.df2 = group_by(matching[,c("cluster_f", "graphst.f")], cluster_f, graphst.f, .drop=FALSE) %>% 
  tally() %>% mutate(log10.n=log10(n+1))

ggplot(plot.df2, aes(x=graphst.f, y=cluster_f, fill=log10.n))+
  geom_tile()+
  scale_fill_gradient2(low="white", mid="skyblue",high="black", midpoint=2)+#median(plot.df2$log10.n))+
  #scale_fill_viridis_c(option="F")#, limits=c(0,1), breaks=c(0,.25,.5,.75,1))+
  labs(fill="log10(# spots)", y="PRECAST cluster", x="GraphST cluster")+
  theme_minimal()+theme(text=element_text(size=16), legend.title=element_text(size=12),
                        legend.text=element_text(size=10))


################## Manual annotation

importCSV.names = grep("csv$", list.files("processed-data/manual_annotation_csv/"), value=T)[1:15]
importCSV = paste0("processed-data/manual_annotation_csv/",importCSV.names)
names(importCSV) = importCSV.names
test = do.call(rbind, lapply(importCSV, read.csv)) %>% tibble::rownames_to_column(var="source")
setdiff(unique(paste(spe$slide, spe$array, sep="_")), unique(test$sample_id))
#ok, even if I recompile the results this slide is still missing


load("processed-data/manual_annotation_csv/compiled_annotation_before_match.Rdata")
colSums(is.na(csv2))
table(csv2$sample_id)
setdiff(unique(paste(spe$slide, spe$array, sep="_")), unique(csv2$sample_id))
#"V12F14-051_C1" "V12F14-051_D1" "V12F14-051_A1" "V12F14-051_B1"
length(setdiff(csv2$spot_name, spe$key))


matching = left_join(as.data.frame(colData(spe)[,c(1:4,15:18,48:51,54)]), csv2[,2:3], by=c("key"="spot_name"))
colSums(is.na(matching))
table(as.character(matching[is.na(matching$ManualAnnotation),"sample_id"]))
match2 = filter(matching, !is.na(ManualAnnotation)) %>% 
  mutate(ManualAnnotation.f=factor(ManualAnnotation, levels=c("THAL","CTX","SUB","PCL-CA1","PCL-CA3","CA4","GCL",
                                                              "ML","SR","SL","SGZ","SLM","SO",
                                                              "WM","CP")))


plot.df = group_by(match2[,c("cluster_f", "ManualAnnotation.f")], cluster_f) %>% 
  add_tally(name="total.precast") %>%
  group_by(cluster_f, ManualAnnotation.f, total.precast, .drop=FALSE) %>% tally() %>%
  mutate(prop.precast=n/total.precast, prop.precast=ifelse(is.na(prop.precast), 0, prop.precast))


ggplot(plot.df, aes(x=ManualAnnotation.f, y=cluster_f, fill=prop.precast))+
  geom_tile()+
  scale_fill_viridis_c(option="F", limits=c(0,1), breaks=c(0,.25,.5,.75,1))+
  labs(fill="prop. PRECAST\ndomain nuclei", y="PRECAST domain", x="Manual cluster")+
  theme_minimal()+theme(text=element_text(size=16), legend.title=element_text(size=12),
                        legend.text=element_text(size=10))

plot.df2 = group_by(match2[,c("cluster_f", "ManualAnnotation.f")], cluster_f, ManualAnnotation.f, .drop=FALSE) %>% 
  tally() %>% mutate(log10.n=log10(n+1))

ggplot(plot.df2, aes(x=ManualAnnotation.f, y=cluster_f, fill=log10.n))+
  geom_tile()+
  scale_fill_gradient2(low="white", mid="skyblue",high="black", midpoint=2)+#median(plot.df2$log10.n))+
  #scale_fill_viridis_c(option="F")#, limits=c(0,1), breaks=c(0,.25,.5,.75,1))+
  labs(fill="log10(# spots)", y="PRECAST cluster", x="Manual annotation")+
  theme_minimal()+theme(text=element_text(size=16), legend.title=element_text(size=12),
                        legend.text=element_text(size=10))#,
                        axis.text.x=element_text(angle=90, hjust=1))
