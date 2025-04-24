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


#alt SLM.SGZ assignment if used k=16
filter(as.data.frame(colData(spe)), domain=="SLM.SGZ") %>% group_by(PRECAST_k18) %>% tally()
#7820 spots total
filter(as.data.frame(colData(spe)), domain=="SLM.SGZ") %>% group_by(PRECAST_k16) %>% tally()
#4825 spots in cluster 5 (SR.SLM)
filter(as.data.frame(colData(spe)), domain=="SLM.SGZ") %>% group_by(PRECAST_k16) %>% tally() %>%
  mutate(perc=n/7820)
table(spe$PRECAST_k16)
#cluster 5: 37003 spots total

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



hm1 = plotGroupedHeatmap(spe, features=genes, group="k16.f", swap_rownames="gene_name",
                   center=T, scale=T, zlim=c(-2,4),
                   cluster_rows=F, cluster_cols=F, border_col="white",
                   color=RColorBrewer::brewer.pal(9,"Greys"),
                   angle_col="0", fontsize=10, fontsize_row=12, silent=T)
hm2 = plotGroupedHeatmap(spe, features=genes, group="k17.f", swap_rownames="gene_name",
                         center=T, scale=T, zlim=c(-2,4),
                         cluster_rows=F, cluster_cols=F, border_col="white",
                         color=RColorBrewer::brewer.pal(9,"Greys"),
                         angle_col="0", fontsize=10, fontsize_row=12, silent=T)
hm3 = plotGroupedHeatmap(spe, features=genes, group="k18.f", swap_rownames="gene_name",
                         center=T, scale=T, zlim=c(-2,4),
                         cluster_rows=F, cluster_cols=F, border_col="white",
                         color=RColorBrewer::brewer.pal(9,"Greys"),
                         angle_col="0", fontsize=10, fontsize_row=12, silent=T)
ggsave(file="plots/revision/supp_PRECAST-compare-k_heatmap.pdf", gridExtra::grid.arrange(hm1[[4]], hm2[[4]], hm3[[4]], ncol=3),
       height=4, width=12, units="in")

################ make pies
#for pies need manual annotation
load("processed-data/manual_annotation_csv/compiled_annotation_before_match.Rdata")
#one brain donor is missing
table(csv2$sample_id)
setdiff(unique(paste(spe$slide, spe$array, sep="_")), unique(csv2$sample_id))
#"V12F14-051_C1" "V12F14-051_D1" "V12F14-051_A1" "V12F14-051_B1"
length(setdiff(csv2$spot_name, spe$key))

#add donor
test = read.csv("processed-data/manual_annotation_csv/spatialLIBD_ManualAnnotation_2023-04-12_Br2720_all.csv") %>%
  mutate(sample_id= as.character(factor(sample_id, levels=c("Br2720_A1","Br2720_B1","Br2720_C1","Br2720_D1"),
                                        labels=c("V12F14-051_A1","V12F14-051_B1","V12F14-051_C1","V12F14-051_D1"))),
         spot_name=paste(spot_name, sample_id, sep="_"))

#combine
csv2 = rbind(csv2, test)

#a handful of un-annotated manual spots so need to subset
dim(spe) #31483 150917
spe$spot_id = paste(sapply(strsplit(rownames(colData(spe)), "_"), function(x) x[[1]]), spe$slide, spe$array, sep="_")
spe.manual = spe[,spe$spot_id %in% csv2$spot_name]
dim(spe.manual) #31483 150520

rownames(csv2) = csv2$spot_name
csv2 = csv2[spe.manual$spot_id,]
identical(spe.manual$spot_id, csv2$spot_name)

spe.manual$ManualAnnotation = csv2$ManualAnnotation

man.pal = c("THAL"="#1e1eff","CTX"="#5ffffb", "SUB"="#add294", "PCL-CA1"="#00dc00", 
             "PCL-CA3"="#00a000", "CA4"="#B0BF1A",
             "GCL"="#005000", "SGZ"="#dfa56e", "ML"="#c1c1c1", "SL"="#444444", "SR"="#828E84", "SLM"="tan4",
             "SO"="#A698AE", "WM"="#ff3ffc", "CP"="#00006a")

p1 <- ggplot(as.data.frame(colData(spe.manual)), aes(x="", fill=ManualAnnotation))+
  geom_bar(stat="count", position="fill", color="white", linewidth=.3)+
  scale_fill_manual(values=man.pal)+
  facet_wrap(vars(k16.f), ncol=5)+theme_void()+
  coord_polar(theta = "y")+labs(title="PRECAST k=16 clusters")+
  theme(aspect.ratio=1, legend.position="none", strip.text=element_text(size=12))
p2 <- ggplot(as.data.frame(colData(spe.manual)), aes(x="", fill=ManualAnnotation))+
  geom_bar(stat="count", position="fill", color="white", linewidth=.3)+
  scale_fill_manual(values=man.pal)+
  facet_wrap(vars(k17.f), ncol=5)+theme_void()+
  coord_polar(theta = "y")+labs(title="PRECAST k=17 clusters")+
  theme(aspect.ratio=1, legend.position="none", strip.text=element_text(size=12))
p3 <- ggplot(as.data.frame(colData(spe.manual)), aes(x="", fill=ManualAnnotation))+
  geom_bar(stat="count", position="fill", color="white", linewidth=.3)+
  scale_fill_manual(values=man.pal)+
  facet_wrap(vars(k18.f), ncol=5)+theme_void()+
  coord_polar(theta = "y")+labs(title="PRECAST k=18 cluster")+
  theme(aspect.ratio=1, legend.position="none", strip.text=element_text(size=12))

pdf(file="plots/revision/supp_PRECAST-compare-k_pie.pdf", height=4, width=12)
gridExtra::grid.arrange(p1+theme(panel.spacing = unit(6, "pt"), plot.margin=unit(c(6,12,6,12),"pt")), 
                        p2+theme(panel.spacing = unit(6, "pt"), plot.margin=unit(c(6,12,6,12),"pt")), 
                        p3+theme(panel.spacing = unit(6, "pt"), plot.margin=unit(c(6,12,6,12),"pt")), ncol=3)
dev.off()

#legend for pies
pl <- ggplot(cbind.data.frame("y"=c(rep(1,4), rep(2,3), rep(3,6), rep(4,2)), 
                        "x"=c(seq(from=1, by=.2, length.out=4),
                              #seq(from=1, by=.2, length.out=3),
                              c(1.4,1.0,1.2),#swap PCL-CA3 to end
			      seq(from=1, by=.2, length.out=6),
                              seq(from=1, by=.2, length.out=2)),
                        "labels"=names(man.pal)),
       aes(x=x, y=y, color=factor(labels, levels=names(man.pal))))+
  geom_point(size=2)+scale_color_manual(values=man.pal)+
  geom_text(aes(label=labels), nudge_x = .03, hjust=0, color="black", size=3)+
  scale_y_reverse()+scale_x_continuous(expand=expansion(add=c(.1,.2)))+
  theme_void()+
  theme(legend.position="none",
        plot.margin=unit(c(.2,.2,.2,.2),"cm"))

ggsave("plots/revision/supp_PRECAST-compare-k_pie-legend.pdf", pl, height=1, width=3)

################ spot plots

test = spe[,spe$sample_id == "Br3942_V11L05-333_D1"]
load(file=here::here('plots','spatial_palette_final.rda'))

k16.pal = c("10"="#5ffffb","8"="#99ff99","1"="#61963d", "14"="#add294", "9"="#00dc00",
                "3"="#85FF33", "7"="#00a000", "16"="#B0BF1A", "4"= "#005000", "2"= "#c1c1c1",
                "13"="#444444", "5"="#777777", 
                "6"="#ff3ffc", "15"="#7a007a", "12"="#00006a", "11"="#1e1eff")

k17.k18.pal = c("10"="#5ffffb","8"="#99ff99","1"="#61963d", "14"="#add294", "9"="#00dc00",
            "3"="#85FF33", "7"="#00a000", "17"="#B0BF1A", "4"= "#005000", "2"= "#c1c1c1",
            "13"="#444444", "5"="#777777", "15"="#dfa56e", 
            "6"="#ff3ffc", "16"="#7a007a", "18"="#ff80fe", "12"="#00006a", "11"="#1e1eff")

p1 <- plotSpots(test, annotate="k16.f", sample_id="sample_id",
          point_size=.5)+
  scale_color_manual(values=k16.pal)+
  guides(colour = guide_legend(ncol = 2))+
  theme_void()+labs(color="k=16")+
  theme(strip.text=element_blank(), legend.key.size = unit(12, "pt"),
        legend.title=element_text(size=12), legend.text=element_text(size=10),
        legend.box.spacing= unit(0, "pt"), legend.margin=margin(0,2,0,0, "pt"))

p2 <- plotSpots(test, annotate="k17.f", sample_id="sample_id",
                point_size=.5)+
  scale_color_manual(values=k17.k18.pal)+
  guides(colour = guide_legend(ncol = 2))+
  theme_void()+labs(color="k=17")+
  theme(strip.text=element_blank(), legend.key.size = unit(12, "pt"),
        legend.title=element_text(size=12), legend.text=element_text(size=10),
        legend.box.spacing= unit(0, "pt"), legend.margin=margin(0,2,0,0, "pt"))

p3 <- plotSpots(test, annotate="k18.f", sample_id="sample_id",
                point_size=.5)+
  scale_color_manual(values=k17.k18.pal)+
  guides(colour = guide_legend(ncol = 2))+
  theme_void()+labs(color="k=18")+
  theme(strip.text=element_blank(), legend.key.size = unit(12, "pt"),
        legend.title=element_text(size=12), legend.text=element_text(size=10),
        legend.box.spacing= unit(0, "pt"), legend.margin=margin(0,2,0,0, "pt"))

pdf(file="plots/revision/supp_PRECAST-compare-k_spot-plots.pdf", height=4, width=12)
gridExtra::grid.arrange(p1, p2, p3, ncol=3)
dev.off()
