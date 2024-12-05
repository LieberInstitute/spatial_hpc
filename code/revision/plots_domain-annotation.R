library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(scater)
set.seed(123)

load(here::here('processed-data','06_clustering','PRECAST','spe_precast_HE_domain.rda'))

############### CA markers

plotGroupedHeatmap(spe, features=c("MPPED1","FIBCD1","CLMP","RGS14","TSPAN18","NECTIN3","AMPH","CARTPT"), 
                   group="cluster", swap_rownames = "gene_name",
                   scale=T, center=T, cluster_rows=F, fontsize=14, angle_col=90)

############### neuropil markers

nrpl = spe[,spe$broad.domain=="Neuropil"]
rownames(nrpl)[1:4]
plotGroupedHeatmap(nrpl, features=c("MPPED1","FIBCD1","TSPAN18","NECTIN3","AMPH",
                                    "PROX1","SEMA5A",
                                    "GFAP","AQP4","APOE","ID4","SOX2"), 
                   group="domain", swap_rownames = "gene_name",
                   scale=T, center=T, cluster_rows=F, fontsize=14, angle_col=0)

############### manual annotation and k18 spot plots

load("processed-data/manual_annotation_csv/compiled_annotation_before_match.Rdata")
test = read.csv("processed-data/manual_annotation_csv/spatialLIBD_ManualAnnotation_2023-04-12_Br2720_all.csv") %>%
  mutate(sample_id= as.character(factor(sample_id, levels=c("Br2720_A1","Br2720_B1","Br2720_C1","Br2720_D1"),
                                        labels=c("V12F14-051_A1","V12F14-051_B1","V12F14-051_C1","V12F14-051_D1"))),
         spot_name=paste(spot_name, sample_id, sep="_"))
csv2 = rbind(csv2, test)

dim(spe) #31483 150917
spe$spot_id = paste(sapply(strsplit(rownames(colData(spe)), "_"), function(x) x[[1]]), spe$slide, spe$array, sep="_")
spe.manual = spe[,spe$spot_id %in% csv2$spot_name]
dim(spe.manual) #31483 150520

rownames(csv2) = csv2$spot_name
csv2 = csv2[spe.manual$spot_id,]
identical(spe.manual$spot_id, csv2$spot_name)

spe.manual$ManualAnnotation = csv2$ManualAnnotation

sub.manual = spe.manual[,spe.manual$sample_id=="Br3942_V11L05-333_D1"]

man.pal = c("THAL"="#1e1eff","CTX"="#5ffffb", "SUB"="#add294", "PCL-CA1"="#00dc00",
            "PCL-CA3"="#00a000", "CA4"="#B0BF1A",
            "GCL"="#005000", "SGZ"="#dfa56e", "ML"="#c1c1c1", "SL"="#444444", "SR"="#828E84", "SLM"="tan4",
            "SO"="#A698AE", "WM"="#ff3ffc", "CP"="#00006a")

plotSpots(sub.manual, annotate="ManualAnnotation", pal=man.pal, point_size=1.5)+
  theme_void()+theme(legend.position="none", aspect.ratio=.8)

k17.k18.pal = c("10"="#5ffffb","8"="#99ff99","1"="#61963d", "14"="#add294", "9"="#00dc00",
                "3"="#85FF33", "7"="#00a000", "17"="#B0BF1A", "4"= "#005000", "2"= "#c1c1c1",
                "13"="#444444", "5"="#777777", "15"="#dfa56e",
                "6"="#ff3ffc", "16"="#7a007a", "18"="#ff80fe", "12"="#00006a", "11"="#1e1eff")

plotSpots(sub.manual, annotate="PRECAST_k18", pal=k17.k18.pal, point_size=1.5)+
  theme_void()+theme(legend.position="none", aspect.ratio=.8)

