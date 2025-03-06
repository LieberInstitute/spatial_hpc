library(SpatialExperiment)
library(ggspavis)
library(dplyr)
library(ggplot2)

load('processed-data/NMF/spe_nmf_final.rda')

# code compiling RCTD results based on here:
# https://github.com/LieberInstitute/spatial_hpc/blob/main/code/revision_maddy/nmfVsRCTD_HE.R

flist = list.files("processed-data/spot_deconvo/RCTD/2ndRun_newClass_RCTDmarkers/layer", pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE)[c(1:32,41:44)]
r.df = do.call(rbind, lapply(flist, function(x) read.csv(x, row.names=NULL)))


### maddy did this but shouldn't have because it doesn't match spe key and actually removes these samples from the merged dataset
#r.df$key <- gsub("Br2720", "V12F14-051", r.df$X)
#colData(spe)[grep("Br2720",spe$key),]
#dim(r.df) #150917
#dim(spe) #150917
#cdata = merge(colData(spe), r.df[,c("key","Oligo","Astro")])
#dim(cdata) #135640

#table(colData(spe)[setdiff(spe$key, r.df$key),c("brnum","sample_id")])
#colData(spe)[grep("Br2720",spe$key),]

##MADDY'S CORR FROM INCOMPLETE DATASET
#cor(cdata$nmf81, cdata$Astro, use = "complete.obs") 
## 0.8577056
#cor(cdata$nmf44, cdata$Oligo, use = "complete.obs") 
## 0.620673

# JT FIX MERGE
cdata = merge(colData(spe), r.df[,c("X","Oligo","Astro")], by.x="key", by.y="X")
dim(cdata) #150917

cor(cdata$nmf81, cdata$Astro, use = "complete.obs") 
# 0.8609265
cor(cdata$nmf44, cdata$Oligo, use = "complete.obs") 
# 0.6141176

#add to spe colData
rownames(cdata) = cdata$key
identical(rownames(cdata), rownames(colData(spe)))
identical(rownames(cdata[rownames(colData(spe)),]), rownames(colData(spe)))
colData(spe) <- cdata[rownames(colData(spe)),]

#pick Br3942 for examples
load("plots/spatial_palette_final.rda")
sub = spe[,spe$brnum=="Br3942"]
#also make an inset version of single panel
length(unique(sub$array_row)) #78 cols
length(unique(sub$array_col)) #128 rows
# array aspect ratio
128/78 # 1.64

##### oligodendrocytes

t1 <- plotVisium(sub, image=F, spots=T, annotate="Oligo", highlight="domain",
           facets="sample_id", point_size=.8)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black")+
  theme_void()+theme(strip.text=element_blank())
ggsave("plots/revision/Figure4_oligo.png", t1, bg="white", height=5.5, width=5.5, units="in")

t2 <- plotVisium(sub, image=F, spots=T, annotate="nmf44", highlight="domain",
                 facets="sample_id", point_size=.8)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black")+
  theme_void()+theme(strip.text=element_blank())
ggsave("plots/revision/Figure4_nmf44.png", t2, bg="white", height=5.5, width=5.5, units="in")

##inset for oligo
sub_oligo = sub[,sub$array_row<42 & sub$sample_id=="V11L05-333_B1"]
length(unique(sub_oligo$array_row)) #42
length(unique(sub_oligo$array_row))*1.64 # 68-69
sub_oligo = sub_oligo[,sub_oligo$array_col> (128-69)]
length(unique(sub_oligo$array_col)) #68

t1.1 <- plotVisium(sub_oligo, image=F, spots=T, annotate="Oligo", highlight="domain",
                 facets="sample_id", point_size=2)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black", limits=c(min(sub$Oligo), max(sub$Oligo)))+
  theme_void()+theme(strip.text=element_blank(), legend.position="none")
ggsave("plots/revision/Figure4_inset-no-legend_oligo.png", t1.1, bg="white", height=2.5, width=2.5, units="in")

t1.2 <- plotVisium(sub_oligo, image=F, spots=T, annotate="Oligo", highlight="domain",
                   facets="sample_id", point_size=1.5)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black", limits=c(min(sub$Oligo), max(sub$Oligo)))+
  theme_void()+theme(strip.text=element_blank())
ggsave("plots/revision/Figure4_inset_oligo.png", t1.2, bg="white", height=2.5, width=2.5, units="in")

t2.1 <- plotVisium(sub_oligo, image=F, spots=T, annotate="nmf44", highlight="domain",
                   facets="sample_id", point_size=2)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black", limits=c(min(sub$nmf44), max(sub$nmf44)))+
  theme_void()+theme(strip.text=element_blank(), legend.position="none")
ggsave("plots/revision/Figure4_inset-no-legend_nmf44.png", t2.1, bg="white", height=2.5, width=2.5, units="in")

t2.2 <- plotVisium(sub_oligo, image=F, spots=T, annotate="nmf44", highlight="domain",
                   facets="sample_id", point_size=1.5)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black", limits=c(min(sub$nmf44), max(sub$nmf44)))+
  theme_void()+theme(strip.text=element_blank())
ggsave("plots/revision/Figure4_inset_nmf44.png", t2.2, bg="white", height=2.5, width=2.5, units="in")

##### astrocytes

t3 <- plotVisium(sub, image=F, spots=T, annotate="Astro", highlight="domain",
                 facets="sample_id", point_size=.8)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black")+
  theme_void()+theme(strip.text=element_blank())
ggsave("plots/revision/Figure4_astro.png", t3, bg="white", height=5.5, width=5.5, units="in")

t4 <- plotVisium(sub, image=F, spots=T, annotate="nmf81", highlight="domain",
                 facets="sample_id", point_size=.8)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black")+
  theme_void()+theme(strip.text=element_blank())
ggsave("plots/revision/Figure4_nmf81.png", t4, bg="white", height=5.5, width=5.5, units="in")

##inset for astro
sub_astro = sub[,sub$array_row>35 & sub$sample_id=="V11L05-333_B1"]
length(unique(sub_astro$array_row)) #42
length(unique(sub_astro$array_row))*1.64 # 68-69
sub_astro = sub_astro[,sub_astro$array_col<68]
length(unique(sub_astro$array_col)) #68

t3.1 <- plotVisium(sub_astro, image=F, spots=T, annotate="Astro", highlight="domain",
                 facets="sample_id", point_size=2)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black", limits=c(min(sub$Astro), max(sub$Astro)))+
  theme_void()+theme(strip.text=element_blank(), legend.position="none")
ggsave("plots/revision/Figure4_inset-no-legend_astro.png", t3.1, bg="white", height=2.5, width=2.5, units="in")

t3.2 <- plotVisium(sub_astro, image=F, spots=T, annotate="Astro", highlight="domain",
                   facets="sample_id", point_size=1.5)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black", limits=c(min(sub$Astro), max(sub$Astro)))+
  theme_void()+theme(strip.text=element_blank())
ggsave("plots/revision/Figure4_inset_astro.png", t3.2, bg="white", height=2.5, width=2.5, units="in")

t4.1 <- plotVisium(sub_astro, image=F, spots=T, annotate="nmf81", highlight="domain",
                   facets="sample_id", point_size=2)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black", limits=c(min(sub$nmf81), max(sub$nmf81)))+
  theme_void()+theme(strip.text=element_blank(), legend.position="none")
ggsave("plots/revision/Figure4_inset-no-legend_nmf81.png", t4.1, bg="white", height=2.5, width=2.5, units="in")

t4.2 <- plotVisium(sub_astro, image=F, spots=T, annotate="nmf81", highlight="domain",
                   facets="sample_id", point_size=1.5)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black", limits=c(min(sub$nmf81), max(sub$nmf81)))+
  theme_void()+theme(strip.text=element_blank())
ggsave("plots/revision/Figure4_inset_nmf81.png", t4.2, bg="white", height=2.5, width=2.5, units="in")


##### neuronal

t5 <- plotVisium(sub, image=F, spots=T, annotate="nmf13", highlight="domain",
                 facets="sample_id", point_size=.8)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black")+
  theme_void()+theme(strip.text=element_blank())
ggsave("plots/revision/Figure4_nmf13.png", t5, bg="white", height=5.5, width=5.5, units="in")

t6 <- plotVisium(sub, image=F, spots=T, annotate="nmf7", highlight="domain",
                 facets="sample_id", point_size=.8)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black")+
  theme_void()+theme(strip.text=element_blank())
ggsave("plots/revision/Figure4_nmf7.png", t6, bg="white", height=5.5, width=5.5, units="in")


##inset for neuronal
sub_nrn = sub[,sub$array_row>35 & sub$sample_id=="V11L05-333_A1"]
length(unique(sub_nrn$array_row)) #42
length(unique(sub_nrn$array_row))*1.64 # 68-69
sub_nrn = sub_nrn[,sub_nrn$array_col<68]
length(unique(sub_nrn$array_col)) #68

t5.1 <- plotVisium(sub_nrn, image=F, spots=T, annotate="nmf13", highlight="domain",
                 facets="sample_id", point_size=2)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black", limits=c(min(sub$nmf13), max(sub$nmf13)))+
  theme_void()+theme(strip.text=element_blank(), legend.position="none")
ggsave("plots/revision/Figure4_inset-no-legend_nmf13.png", t5.1, bg="white", height=2.5, width=2.5, units="in")

t5.2 <- plotVisium(sub_nrn, image=F, spots=T, annotate="nmf13", highlight="domain",
                   facets="sample_id", point_size=1.5)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black", limits=c(min(sub$nmf13), max(sub$nmf13)))+
  theme_void()+theme(strip.text=element_blank())
ggsave("plots/revision/Figure4_inset_nmf13.png", t5.2, bg="white", height=2.5, width=2.5, units="in")

t6.1 <- plotVisium(sub_nrn, image=F, spots=T, annotate="nmf7", highlight="domain",
                   facets="sample_id", point_size=2)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black", limits=c(min(sub$nmf7), max(sub$nmf7)))+
  theme_void()+theme(strip.text=element_blank(), legend.position="none")
ggsave("plots/revision/Figure4_inset-no-legend_nmf7.png", t6.1, bg="white", height=2.5, width=2.5, units="in")

t6.2 <- plotVisium(sub_nrn, image=F, spots=T, annotate="nmf7", highlight="domain",
                   facets="sample_id", point_size=1.5)+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90", high="black", limits=c(min(sub$nmf7), max(sub$nmf7)))+
  theme_void()+theme(strip.text=element_blank())
ggsave("plots/revision/Figure4_inset_nmf7.png", t6.2, bg="white", height=2.5, width=2.5, units="in")


