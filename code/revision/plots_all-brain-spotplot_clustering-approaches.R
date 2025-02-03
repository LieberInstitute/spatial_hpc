library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(scater)
library(gridExtra)
set.seed(123)

load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))
load("plots/spatial_palette_final.rda")

#update sample_id order
spe$sample_id = factor(spe$sample_id, levels=c("V10B01-086_D1","V10B01-086_C1","V11U08-081_C1","V11U08-081_D1",
                                               "V11L05-333_A1","V11L05-333_B1","V11L05-333_C1","V11L05-333_D1",
                                               "V10B01-085_B1","V10B01-085_A1","V10B01-085_D1","V10B01-085_C1",
                                               "V10B01-086_A1","V10B01-086_B1",
                                               #"V11L05-335_C1","V11L05-335_B1","V11L05-335_A1",
                                               "V11L05-335_A1","V11L05-335_C1","V11L05-335_B1",
                                               "V11U08-084_A1","V11U08-084_B1","V11U08-084_C1","V11U08-084_D1",
                                               "V11A20-297_C1","V11A20-297_D1","V11A20-297_A1","V11L05-335_D1","V11A20-297_B1",
                                               "V11U08-081_A1","V11U08-081_B1",
                                               "V11L05-336_A1","V11L05-336_B1","V11L05-336_C1","V11L05-336_D1",
                                               "V12F14-051_C1","V12F14-051_D1","V12F14-051_A1","V12F14-051_B1"))
spe$facet_column = as.character(factor(spe$sample_id, levels=levels(spe$sample_id),
                                       labels=c("c1","c2","c1","c2",
                                                "c1","c2","c1","c2",
                                                "c1","c2","c1","c2",
                                                "c1","c1",
                                                #"c1","c2","c1",
                                                "c1","c2","c3",
                                                "c1","c2","c1","c2",
                                                "c1","c2","c1","c2","c2",
                                                "c1","c1",
                                                "c1","c2","c1","c2",
                                                "c1","c2","c1","c2")))
spe$facet_row = as.character(factor(spe$sample_id, levels=levels(spe$sample_id),
                                    labels=c("r1","r1","r2","r2",
                                             "r1","r1","r2","r2",
                                             "r1","r1","r2","r2",
                                             "r1","r2",
                                             #"r1","r1","r2",
                                             "r1","r1","r1",
                                             "r1","r1","r2","r2",
                                             "r1","r1","r2","r2","r3",
                                             "r1","r2",
                                             "r1","r1","r2","r2",
                                             "r1","r1","r2","r2")))


#brains <- unique(spe$brnum)
brains <- c("Br3942","Br6522","Br8667","Br2743","Br6423","Br2720","Br6432","Br8492","Br6471","Br8325")
p<-list()
for (j in seq_along(brains)) {
  speb <- spe[, (colData(spe)$brnum == brains[j])]
  speb$sample_id <- droplevels(speb$sample_id)
  p[[j]]<-plotVisium(speb, spots=T, image=F, highlight="domain", annotate="domain",
                     facets="sample_id",point_size=.5)+ 
    facet_grid(rows=vars(facet_row), cols=vars(facet_column))+
    ggtitle(brains[j])+scale_y_continuous(expand=expansion(add=4))+scale_x_continuous(expand=expansion(add=4))+
    scale_color_manual(values=spatial.palette, guide="none")+
    scale_fill_manual(values=spatial.palette, guide="none")+guides(color=NULL, fill=NULL)+
    #theme_bw()+
    theme(#legend.text = element_text(size = 8),
      plot.title = element_text(size = 16, hjust=.5, margin=unit(c(5,0,0,0), "pt")),
      strip.text=element_blank(),
      panel.spacing=unit(5,'points'),plot.margin=unit(c(0,5,0,0), "pt"),
    )
}

lplot = ggplot(cbind.data.frame("y"=seq(spatial.palette), "labels"=names(spatial.palette)),
               aes(x=1, y=y, color=labels))+
  geom_point(size=3)+scale_color_manual(values=spatial.palette)+
  geom_text(aes(label=labels), nudge_x = .06, hjust=0, color="black", size=4)+
  scale_y_reverse()+scale_x_continuous(expand=expansion(add=c(.1,.4)))+
  theme_void()+theme(legend.position="none",
                     plot.margin=unit(c(2,.5,2,.5),"cm"))

laymat=rbind(c(1,1,2,2,3,3),c(1,1,2,2,3,3),
             c(4,4,5,5,6,6),c(4,4,5,5,6,6),
             c(7,NA,8,10,10,11),c(7,NA,8,10,10,11),
             c(9,9,9,10,10,11))

pdf(file = "plots/revision/all-brain_spotplot_domain.pdf",
    width=9, height=12)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],lplot, 
             layout_matrix=laymat)
dev.off()

########################## BayesSpace
load("processed-data/06_clustering/BayesSpace/revision/spe_bayes_k18-kmeans-10k.Rdata")

spe$spatial.cluster.f = factor(spe$spatial.cluster, levels=c(16,5,18,3,8,15,4,13,
                                                             14,9,11,10,1,
                                                             7,12,17,6,2
))

spe$sample_id = factor(paste(spe$slide, spe$array, sep="_"), levels=c("V10B01-086_D1","V10B01-086_C1","V11U08-081_C1","V11U08-081_D1",
                                                                      "V11L05-333_A1","V11L05-333_B1","V11L05-333_C1","V11L05-333_D1",
                                                                      "V10B01-085_B1","V10B01-085_A1","V10B01-085_D1","V10B01-085_C1",
                                                                      "V10B01-086_A1","V10B01-086_B1",
                                                                      #"V11L05-335_C1","V11L05-335_B1","V11L05-335_A1",
                                                                      "V11L05-335_A1","V11L05-335_C1","V11L05-335_B1",
                                                                      "V11U08-084_A1","V11U08-084_B1","V11U08-084_C1","V11U08-084_D1",
                                                                      "V11A20-297_C1","V11A20-297_D1","V11A20-297_A1","V11L05-335_D1","V11A20-297_B1",
                                                                      "V11U08-081_A1","V11U08-081_B1",
                                                                      "V11L05-336_A1","V11L05-336_B1","V11L05-336_C1","V11L05-336_D1",
                                                                      "V12F14-051_C1","V12F14-051_D1","V12F14-051_A1","V12F14-051_B1"))
spe$facet_column = as.character(factor(spe$sample_id, levels=levels(spe$sample_id),
                                       labels=c("c1","c2","c1","c2",
                                                "c1","c2","c1","c2",
                                                "c1","c2","c1","c2",
                                                "c1","c1",
                                                #"c1","c2","c1",
                                                "c1","c2","c3",
                                                "c1","c2","c1","c2",
                                                "c1","c2","c1","c2","c2",
                                                "c1","c1",
                                                "c1","c2","c1","c2",
                                                "c1","c2","c1","c2")))
spe$facet_row = as.character(factor(spe$sample_id, levels=levels(spe$sample_id),
                                    labels=c("r1","r1","r2","r2",
                                             "r1","r1","r2","r2",
                                             "r1","r1","r2","r2",
                                             "r1","r2",
                                             #"r1","r1","r2",
                                             "r1","r1","r1",
                                             "r1","r1","r2","r2",
                                             "r1","r1","r2","r2","r3",
                                             "r1","r2",
                                             "r1","r1","r2","r2",
                                             "r1","r1","r2","r2")))

table(colData(spe)[,c("domain","spatial.cluster.f")])

#bayes.palette = c("13"="#005000", "15"="#00a000", "4"="#00a000", "8"="#00dc00",#DG and CA
#                  "3"="#add294", "18"="#61963d", "5"="#99ff99", "16"="#5ffffb",#Sub, RHP, gaba
#                  "10"="#444444", "14"="#c1c1c1", "9"="#777777", "11"="#777777", "1"="#dfa56e",#neuropil
#                  "7"="#ff3ffc", "12"="#7a007a", "17"="#ff80fe", "6"="#1e1eff", "2"="#00006a")
bayes.palette = c("13"="#005000", "15"="#B0BF1A", ################### I CHANGED THIS BECAUSE IT BETTER REFLECTS THE IDENTITY
                  "4"="#00a000", "8"="#00dc00",
                  "3"="#add294", "18"="#61963d", "5"="#99ff99", "16"="#5ffffb",#Sub, RHP, gaba
                  "10"="#444444", "9"="#777777",
                  "11"="#708090", ################### I CHANGED THIS COLOR TO BE SEPARATE FROM 9 WHICH IT WAS ORIGINALLY AND I SHOULD UPDATE IN ANY NECESSARY PLOT 
                  "14"="#c1c1c1",
                  "1"="#dfa56e",#neuropil
                  "7"="#ff3ffc", "12"="#7a007a", "17"="#ff80fe", "6"="#1e1eff", "2"="#00006a")
setdiff(names(bayes.palette), unique(spe$spatial.cluster.f))
setdiff(unique(spe$spatial.cluster.f),names(bayes.palette))


brains <- c("Br3942","Br6522","Br8667","Br2743","Br6423","Br2720","Br6432","Br8492","Br6471","Br8325")
p<-list()
for (j in seq_along(brains)) {
  speb <- spe[, (colData(spe)$brnum == brains[j])]
  speb$sample_id <- droplevels(speb$sample_id)
  p[[j]]<-plotVisium(speb, spots=T, image=F, highlight="spatial.cluster.f", annotate="spatial.cluster.f",
                     facets="sample_id",point_size=.5)+ 
    facet_grid(rows=vars(facet_row), cols=vars(facet_column))+
    ggtitle(brains[j])+scale_y_continuous(expand=expansion(add=4))+scale_x_continuous(expand=expansion(add=4))+
    scale_color_manual(values=bayes.palette, guide="none")+
    scale_fill_manual(values=bayes.palette, guide="none")+guides(color=NULL, fill=NULL)+
    #theme_bw()+
    theme(#legend.text = element_text(size = 8),
      plot.title = element_text(size = 16, hjust=.5, margin=unit(c(5,0,0,0), "pt")),
      strip.text=element_blank(),
      panel.spacing=unit(5,'points'),plot.margin=unit(c(0,5,0,0), "pt"),
    )
}

lplot = ggplot(cbind.data.frame("y"=seq(bayes.palette), "labels"=names(bayes.palette)),
               aes(x=1, y=y, color=labels))+
  geom_point(size=3)+scale_color_manual(values=bayes.palette)+
  geom_text(aes(label=labels), nudge_x = .06, hjust=0, color="black", size=4)+
  scale_y_reverse()+scale_x_continuous(expand=expansion(add=c(.1,.4)))+
  theme_void()+theme(legend.position="none",
                     plot.margin=unit(c(2,.5,2,.5),"cm"))

laymat=rbind(c(1,1,2,2,3,3),c(1,1,2,2,3,3),
             c(4,4,5,5,6,6),c(4,4,5,5,6,6),
             c(7,NA,8,10,10,11),c(7,NA,8,10,10,11),
             c(9,9,9,10,10,11))

pdf(file = "plots/revision/all-brain_spotplot_BayesSpace-recolored.pdf",
    width=9, height=12)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],lplot, 
             layout_matrix=laymat)
dev.off()

########################## GraphST
graphst = read.csv("processed-data/06_clustering/GraphST/hpc_allSamples_graphst_k16.csv") %>%
  mutate(spot_id = gsub("1-[0-9]*", "1_", spot_id))
graphst$spot_id = paste0(graphst$spot_id, graphst$captureArea)
colSums(is.na(graphst))

#fix spot codes in colData
spe$spot_id = paste(sapply(strsplit(rownames(colData(spe)), "_"), function(x) x[[1]]), spe$slide, spe$array, sep="_")

#reorder graphst
rownames(graphst) = graphst$spot_id
graphst = graphst[spe$spot_id,]
identical(graphst$spot_id, spe$spot_id)

spe$graphst = as.character(graphst$cluster_lamb_0_1_and_1)

#graphst.palette = c("11"="#005000","2"="#00a000","9"="#00a000","3"="#00dc00",
#                    "15"="#add294","14"="#61963d","7"="#99ff99","4"="#5ffffb",
#                    "8"="#c1c1c1","6"="#777777","5"="#dfa56e",
#                    "13"="#ff3ffc","12"="#7a007a","16"="#ff80fe","10"="#1e1eff","1"="#00006a")
graphst.palette = c("11"="#005000","2"="#B0BF1A", ############### I CHANGED THIS  BECAUSE IT BETTER REFLECTS HIS IDENTITY
                    "9"="#00a000", #CA3
                    "7"="#00A36C",############### I CHANGED THIS COLOR TO BE SEPARATE FROM ANY OTHER GREEN BECAUSE IT IS A MIX OF ALL CA PYRAMIDALS AND I SHOULD UPDATE IN ANY NECESSARY PLOT 
                    "3"="#00dc00",
                    "15"="#add294","14"="#61963d",
                    "4"="#5ffffb",
                    "6"="#777777","8"="#c1c1c1","5"="#dfa56e",
                    "13"="#ff3ffc","12"="#7a007a","16"="#ff80fe","10"="#1e1eff","1"="#00006a")

setdiff(names(graphst.palette), unique(spe$graphst))
setdiff(unique(spe$graphst),names(graphst.palette))

brains <- c("Br3942","Br6522","Br8667","Br2743","Br6423","Br2720","Br6432","Br8492","Br6471","Br8325")
p<-list()
for (j in seq_along(brains)) {
  speb <- spe[, (colData(spe)$brnum == brains[j])]
  speb$sample_id <- droplevels(speb$sample_id)
  p[[j]]<-plotVisium(speb, spots=T, image=F, highlight="graphst", annotate="graphst",
                     facets="sample_id",point_size=.5)+ 
    facet_grid(rows=vars(facet_row), cols=vars(facet_column))+
    ggtitle(brains[j])+scale_y_continuous(expand=expansion(add=4))+scale_x_continuous(expand=expansion(add=4))+
    scale_color_manual(values=graphst.palette, guide="none")+
    scale_fill_manual(values=graphst.palette, guide="none")+guides(color=NULL, fill=NULL)+
    #theme_bw()+
    theme(#legend.text = element_text(size = 8),
      plot.title = element_text(size = 16, hjust=.5, margin=unit(c(5,0,0,0), "pt")),
      strip.text=element_blank(),
      panel.spacing=unit(5,'points'),plot.margin=unit(c(0,5,0,0), "pt"),
    )
}

lplot = ggplot(cbind.data.frame("y"=seq(graphst.palette), "labels"=names(graphst.palette)),
               aes(x=1, y=y, color=labels))+
  geom_point(size=3)+scale_color_manual(values=graphst.palette)+
  geom_text(aes(label=labels), nudge_x = .06, hjust=0, color="black", size=4)+
  scale_y_reverse()+scale_x_continuous(expand=expansion(add=c(.1,.4)))+
  theme_void()+theme(legend.position="none",
                     plot.margin=unit(c(2,.5,2,.5),"cm"))

laymat=rbind(c(1,1,2,2,3,3),c(1,1,2,2,3,3),
             c(4,4,5,5,6,6),c(4,4,5,5,6,6),
             c(7,NA,8,10,10,11),c(7,NA,8,10,10,11),
             c(9,9,9,10,10,11))

pdf(file = "plots/revision/all-brain_spotplot_GraphST-recolored.pdf",
    width=9, height=12)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],lplot, 
             layout_matrix=laymat)
dev.off()

########################## Manual annotation
load("processed-data/manual_annotation_csv/compiled_annotation_before_match.Rdata")
test = read.csv("processed-data/manual_annotation_csv/spatialLIBD_ManualAnnotation_2023-04-12_Br2720_all.csv") %>%
  mutate(sample_id= as.character(factor(sample_id, levels=c("Br2720_A1","Br2720_B1","Br2720_C1","Br2720_D1"),
                                        labels=c("V12F14-051_A1","V12F14-051_B1","V12F14-051_C1","V12F14-051_D1"))),
         spot_name=paste(spot_name, sample_id, sep="_")) 
csv2 = rbind(csv2, test)

dim(spe) #31483 150917
spe.manual = spe[,spe$spot_id %in% csv2$spot_name]
dim(spe.manual) #31483 150520

#reorg test
rownames(csv2) = csv2$spot_name
csv2 = csv2[spe.manual$spot_id,]
identical(spe.manual$spot_id, csv2$spot_name)

spe.manual$ManualAnnotation = csv2$ManualAnnotation

#manual.palette = c("GCL"="#005000", "CA4"="#B0BF1A", "PCL-CA3"="#00a000", "PCL-CA1"="#00dc00",
#                   "SUB"="#add294", "CTX"="#99ff99", "THAL"="#1e1eff", 
#                   "SL"="#444444", "SO"="#A698AE", "SR"="#828E84", "ML"="#c1c1c1", "SLM"="tan4","SGZ"="#dfa56e", 
#                   "WM"="#ff3ffc", "CP"="#00006a")
manual.palette = c("GCL"="#005000", "CA4"="#B0BF1A", "PCL-CA3"="#00a000", "PCL-CA1"="#00dc00",
                   "SUB"="#add294", "CTX"="#99ff99", 
                   #"THAL"="#1e1eff", ## i need to change this because this is vascular color for everyone else 
                   "THAL"="#6495ED", ## changed
                   "SL"="#444444", "SO"="#A698AE", "SR"="#828E84", "SLM"="tan4","ML"="#c1c1c1", "SGZ"="#dfa56e", 
                   "WM"="#ff3ffc", "CP"="#00006a")
setdiff(names(manual.palette), spe.manual$ManualAnnotation)
setdiff(spe.manual$ManualAnnotation,names(manual.palette))

brains <- c("Br3942","Br6522","Br8667","Br2743","Br6423","Br2720","Br6432","Br8492","Br6471","Br8325")
p<-list()
for (j in seq_along(brains)) {
  speb <- spe.manual[, (colData(spe.manual)$brnum == brains[j])]
  speb$sample_id <- droplevels(speb$sample_id)
  p[[j]]<-plotVisium(speb, spots=T, image=F, highlight="ManualAnnotation", annotate="ManualAnnotation",
                     facets="sample_id",point_size=.5)+ 
    facet_grid(rows=vars(facet_row), cols=vars(facet_column))+
    ggtitle(brains[j])+scale_y_continuous(expand=expansion(add=4))+scale_x_continuous(expand=expansion(add=4))+
    scale_color_manual(values=manual.palette, guide="none")+
    scale_fill_manual(values=manual.palette, guide="none")+guides(color=NULL, fill=NULL)+
    #theme_bw()+
    theme(#legend.text = element_text(size = 8),
      plot.title = element_text(size = 16, hjust=.5, margin=unit(c(5,0,0,0), "pt")),
      strip.text=element_blank(),
      panel.spacing=unit(5,'points'),plot.margin=unit(c(0,5,0,0), "pt"),
    )
}

lplot = ggplot(cbind.data.frame("y"=seq(manual.palette), "labels"=names(manual.palette)),
               aes(x=1, y=y, color=labels))+
  geom_point(size=3)+scale_color_manual(values=manual.palette)+
  geom_text(aes(label=labels), nudge_x = .06, hjust=0, color="black", size=4)+
  scale_y_reverse()+scale_x_continuous(expand=expansion(add=c(.1,.4)))+
  theme_void()+theme(legend.position="none",
                     plot.margin=unit(c(2,.5,2,.5),"cm"))

laymat=rbind(c(1,1,2,2,3,3),c(1,1,2,2,3,3),
             c(4,4,5,5,6,6),c(4,4,5,5,6,6),
             c(7,NA,8,10,10,11),c(7,NA,8,10,10,11),
             c(9,9,9,10,10,11))

pdf(file = "plots/revision/all-brain_spotplot_ManualAnnotation-recolored.pdf",
    width=9, height=12)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],lplot, 
             layout_matrix=laymat)
dev.off()
