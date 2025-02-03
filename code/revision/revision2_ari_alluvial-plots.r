library(SpatialExperiment)
library(ggspavis)
library(dplyr)
library(ggplot2)
library(ggsankey)
library(mclust)
#merging of clusters based on code from https://github.com/LieberInstitute/spatial_hpc/blob/main/code/revision/plots_spatial-clustering-approaches.R

#load in spe object with both precast and bayesspace
load("processed-data/06_clustering/BayesSpace/revision/spe_bayes_k18-kmeans-10k.Rdata")

#fix spot id so it matches with graphst
spe$spot_id = paste(sapply(strsplit(rownames(colData(spe)), "_"), function(x) x[[1]]), spe$slide, spe$array, sep="_")
spe$BayesSpace = spe$spatial.cluster
#match in graphst labels
graphst = read.csv("processed-data/06_clustering/GraphST/hpc_allSamples_graphst_k16.csv") %>%
  mutate(spot_id = gsub("1-[0-9]*", "1_", spot_id))
graphst$spot_id = paste0(graphst$spot_id, graphst$captureArea)

clus.df = left_join(as.data.frame(colData(spe))[,c("spot_id","brnum","PRECAST_k18","cluster","domain","BayesSpace")] %>%
                       mutate(brnum=as.character(brnum)),
                     graphst, by=c("spot_id", "brnum"="donor_id"))
colSums(is.na(clus.df))
clus.df$GraphST = clus.df$cluster_lamb_0_1_and_1
#add in manual annotation
load("processed-data/manual_annotation_csv/compiled_annotation_before_match.Rdata")
test = read.csv("processed-data/manual_annotation_csv/spatialLIBD_ManualAnnotation_2023-04-12_Br2720_all.csv") %>%
  mutate(sample_id= as.character(factor(sample_id, levels=c("Br2720_A1","Br2720_B1","Br2720_C1","Br2720_D1"),
                                        labels=c("V12F14-051_A1","V12F14-051_B1","V12F14-051_C1","V12F14-051_D1"))),
         spot_name=paste(spot_name, sample_id, sep="_")) 
csv2 = rbind(csv2, test)

match_manual = left_join(clus.df[,c("spot_id","PRECAST_k18","domain","cluster","BayesSpace","GraphST")], 
                         csv2[,2:3], by=c("spot_id"="spot_name"))
colSums(is.na(match_manual)) #should only be 397 missing/NA for ManualAnnotation
match_manual = filter(match_manual, !is.na(ManualAnnotation))


#copied from https://github.com/LieberInstitute/spatial_hpc/blob/main/code/revision/plots_PRECAST-k-compare.R
manual.palette = c("GCL"="#005000", "CA4"="#B0BF1A", "PCL-CA3"="#00a000", "PCL-CA1"="#00dc00",
                   "SUB"="#add294", "CTX"="#99ff99", 
                   #"THAL"="#1e1eff", ## i need to change this because this is vascular color for everyone else 
                   "THAL"="#6495ED", ## changed
                   "SL"="#444444", "SO"="#A698AE", "SR"="#828E84", "SLM"="tan4","ML"="#c1c1c1", "SGZ"="#dfa56e", 
                   "WM"="#ff3ffc", "CP"="#00006a")


#pairwise with precast
#k18.pal reordered from k17.18.pal here: https://github.com/LieberInstitute/spatial_hpc/blob/main/code/revision/plots_PRECAST-k-compare.R
k18.pal = c("4"= "#005000",
            "17"="#B0BF1A","7"="#00a000",#two ca3 clusters
            "9"="#00dc00","3"="#85FF33",#two ca1 clusters
            "14"="#add294",#sub
            "1"="#61963d",#sub.rhp
            "8"="#99ff99",#rhp/cortex
            "10"="#5ffffb",#GABA
            "13"="#444444","5"="#777777","2"= "#c1c1c1","15"="#dfa56e", #SL.SR, SR.SLM, ML, SLM.SGZ
            "6"="#ff3ffc", "16"="#7a007a", "18"="#ff80fe", "11"="#1e1eff", "12"="#00006a" 
            )


manual_precast = make_long(match_manual, ManualAnnotation, PRECAST_k18)
manual_precast = mutate(manual_precast, node=factor(node, levels=c(names(manual.palette)[length(manual.palette):1],
                                                                       names(k18.pal)[length(k18.pal):1])
                                                    )
                        )


man.pre.palette = c(manual.palette, k18.pal)

p1 <- ggplot(manual_precast, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=node))+
  geom_alluvial(node.color="white", linewidth=.5)+
  scale_fill_manual(values=man.pre.palette)+
  scale_x_discrete(expand=expansion(add=c(.3,.2)))+
  scale_y_continuous(labels=c("0","50k","100k","150k"))+
  #geom_alluvial_label(aes(label=node), hjust=1, fill="white")+
  labs(title="Manual vs. PRECAST k=18", x="", y="# of spots")+
  theme_bw()+theme(legend.position="none", axis.title.x=element_blank(),
                   axis.text.x=element_text(size=12), axis.title.y=element_text(size=12),
                   panel.grid.major.x = element_blank(), plot.margin=margin(.5,.5,.5,.5, unit="cm"))

#ari
adjustedRandIndex(match_manual$ManualAnnotation, match_manual$PRECAST_k18) #0.1922123
adjustedRandIndex(match_manual$ManualAnnotation, match_manual$domain) #0.2312085
  
#pairwise with GraphST
#palette from https://github.com/LieberInstitute/spatial_hpc/blob/main/code/revision/plots_all-brain-spotplot_clustering-approaches.R
#change order for better organization of sankey plot
graphst.palette = c("11"="#005000","2"="#B0BF1A", ############### I CHANGED THIS  BECAUSE IT BETTER REFLECTS HIS IDENTITY
                    "9"="#00a000", #CA3
                    "7"="#00A36C",############### I CHANGED THIS COLOR TO BE SEPARATE FROM ANY OTHER GREEN BECAUSE IT IS A MIX OF ALL CA PYRAMIDALS AND I SHOULD UPDATE IN ANY NECESSARY PLOT 
                    "3"="#00dc00",
                    "15"="#add294","14"="#61963d",
                    "4"="#5ffffb",
                    "6"="#777777","8"="#c1c1c1","5"="#dfa56e",
                    "13"="#ff3ffc","12"="#7a007a","16"="#ff80fe","10"="#1e1eff","1"="#00006a")

manual_gst = make_long(match_manual, ManualAnnotation, GraphST)
manual_gst = mutate(manual_gst, node=factor(node, levels=c(names(manual.palette)[length(manual.palette):1],
                                                           names(graphst.palette)[length(graphst.palette):1])
                                            )
                    )

man.gst.palette = c(manual.palette, graphst.palette)

p2 <- ggplot(manual_gst, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=node))+
  geom_alluvial(node.color="white", linewidth=.5)+
  scale_fill_manual(values=man.gst.palette)+
  scale_x_discrete(expand=expansion(add=c(.3,.2)))+
  scale_y_continuous(labels=c("0","50k","100k","150k"))+
  #geom_alluvial_label(aes(label=node), hjust=1, fill="white")+
  labs(title="Manual vs. GraphST", x="", y="# of spots")+
  theme_bw()+theme(legend.position="none", axis.title.x=element_blank(),
                   axis.text.x=element_text(size=12), axis.title.y=element_text(size=12),
                   panel.grid.major.x = element_blank(), plot.margin=margin(.5,.5,.5,.5, unit="cm"))

#ari
adjustedRandIndex(match_manual$ManualAnnotation, match_manual$GraphST) #0.11179

#pairwise with BayesSpace
#palette from palette from https://github.com/LieberInstitute/spatial_hpc/blob/main/code/revision/plots_all-brain-spotplot_clustering-approaches.R
#rearrange for better order in sankey plot
bayes.palette = c("13"="#005000", "15"="#B0BF1A", ################### I CHANGED THIS BECAUSE IT BETTER REFLECTS THE IDENTITY
                  "4"="#00a000", "8"="#00dc00",
                  "3"="#add294", "18"="#61963d", "5"="#99ff99", "16"="#5ffffb",#Sub, RHP, gaba
                  "10"="#444444", "9"="#777777",
                  "11"="#708090", ################### I CHANGED THIS COLOR TO BE SEPARATE FROM 9 WHICH IT WAS ORIGINALLY AND I SHOULD UPDATE IN ANY NECESSARY PLOT 
                  "14"="#c1c1c1",
                  "1"="#dfa56e",#neuropil
                  "7"="#ff3ffc", "12"="#7a007a", "17"="#ff80fe", "6"="#1e1eff", "2"="#00006a")


manual_bsp = make_long(match_manual, ManualAnnotation, BayesSpace)
manual_bsp = mutate(manual_bsp, node=factor(node, levels=c(names(manual.palette)[length(manual.palette):1],
                                                           names(bayes.palette)[length(bayes.palette):1])
                                            )
                    )

man.bsp.palette = c(manual.palette, bayes.palette)

p3 <- ggplot(manual_bsp, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=node))+
  geom_alluvial(node.color="white", linewidth=.5)+
  scale_fill_manual(values=man.bsp.palette)+
  scale_x_discrete(expand=expansion(add=c(.3,.2)))+
  scale_y_continuous(labels=c("0","50k","100k","150k"))+
  #geom_alluvial_label(aes(label=node), hjust=1, fill="white")+
  labs(title="Manual vs. BayesSpace", x="", y="# of spots")+
  theme_bw()+theme(legend.position="none", axis.title.x=element_blank(),
                   axis.text.x=element_text(size=12), axis.title.y=element_text(size=12),
                   panel.grid.major.x = element_blank(), plot.margin=margin(.5,.5,.5,.5, unit="cm"))

#ari
adjustedRandIndex(match_manual$ManualAnnotation, match_manual$BayesSpace) #0.1727502
match_manual$bsp_domain= factor(match_manual$BayesSpace, levels=1:18, labels=c(1:3,"CA2.4",5:14,"CA2.4",16:18))
adjustedRandIndex(match_manual$ManualAnnotation, match_manual$bsp_domain) #0.178568


#precast with bayesspace
#renamed k18.pal
tmp = distinct(match_manual, PRECAST_k18, cluster)
rename1 = tmp$cluster
names(rename1) = as.character(tmp$PRECAST_k18)
k18.pal2 = k18.pal 
rename1.order = rename1[names(k18.pal2)]
names(k18.pal2) = rename1.order

prc_bsp = make_long(match_manual, cluster, BayesSpace)
prc_bsp = mutate(prc_bsp, node=factor(node, levels=c(names(k18.pal2)[length(k18.pal2):1],
                                                           names(bayes.palette)[length(bayes.palette):1])
                                      )
                 )

k18.bsp.palette = c(k18.pal2, bayes.palette)

p4 <- ggplot(prc_bsp, aes(x=x, next_x=next_x, node=node, next_node=next_node, fill=node))+
  geom_alluvial(node.color="white", linewidth=.5)+
  scale_fill_manual(values=k18.bsp.palette)+
  scale_x_discrete(expand=expansion(add=c(.3,.2)), labels=c("PRECAST_k18","BayesSpace"))+
  scale_y_continuous(labels=c("0","50k","100k","150k"))+
  #geom_alluvial_label(aes(label=node), hjust=1, fill="white")+
  labs(title="PRECAST k=18 vs. BayesSpace", x="", y="# of spots")+
  theme_bw()+theme(legend.position="none", axis.title.x=element_blank(),
                   axis.text.x=element_text(size=12), axis.title.y=element_text(size=12),
                   panel.grid.major.x = element_blank(), plot.margin=margin(.5,.5,.5,.5, unit="cm"))

#ari
adjustedRandIndex(clus.df$cluster, clus.df$BayesSpace) #0.3983147
clus.df$bsp_domain = factor(clus.df$BayesSpace, levels=1:18, labels=c(1:3,"CA2.4",5:14,"CA2.4",16:18))
adjustedRandIndex(clus.df$domain, clus.df$bsp_domain) #0.4265898

adjustedRandIndex(clus.df$cluster, clus.df$GraphST) #0.2342174

#save alluvial
#dims currently for 1 plot
ggsave(filename="plots/revision/clustering-approaches_alluvial.png", 
       gridExtra::grid.arrange(p1, p2, p3, p4), 
       bg="white", units="in", height=10, width=10)



################# pairwise jaccard
#jaccard coeff
jcoef <- function(x, y) { 
  # x is a named T/F vector for reference cluster ID, with names being spot codes
  # y is a named factor, numeric, or character vector of all comparison cluster IDs, with names being spotcodes
  if(class(y)!="factor") y=as.factor(y)
  x.ids = names(x)[x]
  y.list = levels(y)
  names(y.list) = y.list
  #for all levels of comparison factor, return JC with reference
  sapply(y.list, function(z) {
    z.ids = names(y)[y==z]
    length(intersect(x.ids, z.ids))/length(union(x.ids, z.ids))
  })
}

pairwise_jc <- function(reference_type, compare_type) {
  #reference type and compare type are columns in match_manual
  x_values = unique(match_manual[[reference_type]])
  #for all values of reference cluster compute JC
  output_list <- lapply(x_values, function(X) {
    test_x = match_manual[[reference_type]]==X
    names(test_x) = match_manual$spot_id
    test_y = match_manual[[compare_type]]
    names(test_y) = match_manual$spot_id
    #for reference cluster ID X, compute JC for all comparison clusters
    jc_output <- jcoef(test_x, test_y)
    cbind.data.frame(ref_type=rep(reference_type, length(jc_output)),
                     ref_clus=rep(X, length(jc_output)),
                     comp_type=rep(compare_type, length(jc_output)),
                     comp_clus=names(jc_output),
                     j.coef=as.numeric(jc_output))
  })
  do.call(rbind, output_list)
}

man.jc.df = bind_rows(pairwise_jc("ManualAnnotation","PRECAST_k18"),
                      pairwise_jc("ManualAnnotation","GraphST"),
                      pairwise_jc("ManualAnnotation","BayesSpace")) %>%
  mutate(comp_type= factor(comp_type, levels=c("PRECAST_k18","GraphST","BayesSpace")),
         ref_clus= factor(ref_clus, levels=names(manual.palette)))


ggplot(man.jc.df, aes(x=comp_type, y=j.coef, label=as.character(comp_clus)))+
  geom_text()+facet_wrap(vars(ref_clus), ncol=5)+#ylim(0,1)+
  labs(title="Manual Annotation", x="spatial clustering methods", y="Jaccard coef.")+
  theme(axis.text.x=element_text(angle=45, hjust=1))

#this plot is just submitted as part of the response to reviewer 
