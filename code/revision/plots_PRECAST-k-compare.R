library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(scater)


load(here::here('processed-data','06_clustering','PRECAST','spe_precast_HE_domain.rda'))

################## load in k=16 and k=17 results
lowk = read.csv("processed-data/revision/PRECAST_k16_k17.csv", row.names=1)
identical(rownames(lowk),colnames(spe))
colData(spe)$PRECAST_k16 = lowk$PRECAST_k16
colData(spe)$PRECAST_k17 = lowk$PRECAST_k17

################## heatmap
genes = c(#"SHOX2",
  #"TCF7L2",
  "GAD2","SST",
  #"SLC17A6",
  #"CDH22","OPRM1","CACNG4",
  "RORB","CUX2",#"SATB2","TLE4",
  "NTS","FN1",
  "MPPED1","FIBCD1",
  "TSPAN18","NECTIN3","AMPH",
  #"SFRP2",
  "PROX1","SEMA5A",
  #"SOX2","SLC1A2",
  "MT-ND5","MT-ATP8","GFAP",
  "MOG","MOBP","PLP1","PRLR","TTR","PDGFRB","TPM2"
)

spe$k16.f = factor(spe$PRECAST_k16, levels=c(10,8,1,
                                             14,3,9,
                                             7,16,
                                             4,2,5,13,
                                             6,15,12,11
))
spe$k17.f = factor(spe$PRECAST_k17, levels=c(10,8,1,
                                             14,3,9,7,17,
                                             4,2,
                                             5,13,15,
                                             6,16,12,11
))
spe$k18.f = factor(spe$PRECAST_k18, levels=c(10,8,1,
                                             14,3,9,7,17,
                                             4,2,
                                             5,13,15,
                                             6,16,18,
                                             12,11
))



plotGroupedHeatmap(spe, features=genes, group="k16.f", swap_rownames="gene_name",
                   center=T, scale=T, zlim=c(-2,4),
                   cluster_rows=F, cluster_cols=F, border_col="white",
                   color=RColorBrewer::brewer.pal(9,"Greys"),
                   angle_col="0")



################ make pies
#for pies need manual annotation
load("processed-data/manual_annotation_csv/compiled_annotation_before_match.Rdata")
#one brain donor is missing
table(csv2$sample_id)
setdiff(unique(paste(spe$slide, spe$array, sep="_")), unique(csv2$sample_id))
#"V12F14-051_C1" "V12F14-051_D1" "V12F14-051_A1" "V12F14-051_B1"
length(setdiff(csv2$spot_name, spe$key))


#matching = left_join(as.data.frame(colData(spe)[,c(1:4,60,15,48,58,59,62)]), csv2[,2:3], by=c("key"="spot_name"))
#precast only (below)
matching = left_join(as.data.frame(colData(spe)[,c(1:4,15,48,52:56)]), csv2[,2:3], by=c("key"="spot_name"))
colSums(is.na(matching))
table(as.character(matching[is.na(matching$ManualAnnotation),"sample_id"]))

match2 = filter(matching, !is.na(ManualAnnotation)) %>% 
  mutate(ManualAnnotation.f=factor(ManualAnnotation, levels=c("THAL","CTX","SUB","PCL-CA1","PCL-CA3","CA4","GCL",
                                                              "SGZ","ML","SL","SR","SLM","SO",
                                                              "WM","CP")))

man.pal = c("THAL"="#1e1eff","CTX"="#5ffffb", "SUB"="#add294", "PCL-CA1"="#00dc00", 
             "PCL-CA3"="#00a000", "CA4"="#B0BF1A",
             "GCL"="#005000", "SGZ"="#dfa56e", "ML"="#c1c1c1", "SL"="#444444", "SR"="#828E84", "SLM"="tan4",
             "SO"="#A698AE", "WM"="#ff3ffc", "CP"="#00006a")

ggplot(match2, aes(x="", fill=ManualAnnotation.f))+
  geom_bar(stat="count", position="fill")+
  scale_fill_manual(values=man.pal)+
  facet_wrap(vars(k18.f), ncol=5)+theme_void()+
  coord_polar(theta = "y")+labs(title="PRECAST k=18 cluster")


################ spot plots

test = spe[,spe$sample_id %in% c("Br3942_V11L05-333_C1","Br3942_V11L05-333_D1")]
load(file=here::here('plots','spatial_palette_final.rda'))

k16.pal = c("10"="#5ffffb","8"="#99ff99","1"="#61963d", "14"="#add294", "9"="#00dc00",
                "3"="#85FF33", "7"="#00a000", "16"="#B0BF1A", "4"= "#005000", "2"= "#c1c1c1",
                "13"="#444444", "5"="#777777", 
                "6"="#ff3ffc", "15"="#7a007a", "12"="#00006a", "11"="#1e1eff")

k17.k18.pal = c("10"="#5ffffb","8"="#99ff99","1"="#61963d", "14"="#add294", "9"="#00dc00",
            "3"="#85FF33", "7"="#00a000", "17"="#B0BF1A", "4"= "#005000", "2"= "#c1c1c1",
            "13"="#444444", "5"="#777777", "15"="#dfa56e", 
            "6"="#ff3ffc", "16"="#7a007a", "18"="#ff80fe", "12"="#00006a", "11"="#1e1eff")

#test$is15_k17 = test$PRECAST_k17==15
#test$is15_k18 = test$PRECAST_k18==15
plotSpots(test, annotate="k16.f", sample_id="sample_id",
          point_size=.5)+
  scale_color_manual(values=k16.pal)+theme_void()+ggtitle("k=16")

plotSpots(test, annotate="k17.f", sample_id="sample_id",
          point_size=.5)+
  scale_color_manual(values=k17.k18.pal)+theme_void()+ggtitle("k=17")

plotSpots(test, annotate="k18.f", sample_id="sample_id",
          point_size=.5)+
  scale_color_manual(values=k17.k18.pal)+theme_void()+ggtitle("k=18")
