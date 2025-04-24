library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(scater)
set.seed(123)

load(here::here('processed-data','06_clustering','PRECAST','spe_precast_HE_domain.rda'))
check_spe <- readRDS("processed-data/spot_deconvo/shared_utilities/spe.rds")
identical(colnames(spe), colnames(check_spe))
spe$nuclei_count = check_spe$count

############### CA markers
nrn = spe[,spe$broad.domain=="Neuron"]
hm1 = plotGroupedHeatmap(nrn, features=c("MPPED1","FIBCD1","CLMP","RGS14","TSPAN18","NECTIN3","AMPH","CARTPT"), 
                   group="cluster", swap_rownames = "gene_name",
                   scale=T, center=T, cluster_rows=F, treeheight_col=15, fontsize=7, angle_col=90,
		   fontsize_row=10, fontsize_col=10)

p1 <- ggplot(filter(as.data.frame(colData(spe)), cluster %in% c("CA2.4.1","CA2.4.2","CA1.1","CA1.2")), 
       aes(x=cluster, y=log2(nuclei_count+1), fill=domain))+
  geom_boxplot(outlier.size = .5, linewidth=.5)+scale_fill_manual(values=c("CA1"="#00dc00","CA2.4"="#00a000"))+
  labs(y="VistoSeg nuclei count (log2)")+
  theme_minimal()+theme(legend.position="none", axis.title.x=element_blank(),
                        axis.title.y=element_text(size=10), axis.text=element_text(size=8))

p2 <- ggplot(filter(as.data.frame(colData(spe)), cluster %in% c("CA2.4.1","CA2.4.2","CA1.1","CA1.2")), 
             aes(x=cluster, y=log10(sum), fill=domain))+
  geom_boxplot(outlier.size = .5, linewidth=.5)+scale_fill_manual(values=c("CA1"="#00dc00","CA2.4"="#00a000"))+
  labs(y="Library size (log10)")+
  theme_minimal()+theme(legend.position="none", axis.title.x=element_blank(),
                        axis.title.y=element_text(size=10), axis.text=element_text(size=8))
  
ggsave(file="plots/revision/supp_cluster-annotation_CA1-CA3-plots.pdf", 
       gridExtra::grid.arrange(hm1[[4]], p1, p2, ncol=3),#layout_matrix=cbind(c(1,1),c(1,1),c(2,3))),
       width=7, height=2.5, units="in"
)
############### neuropil markers

nrpl = spe[,spe$broad.domain=="Neuropil"]
rownames(nrpl)[1:4]
hm2 = plotGroupedHeatmap(nrpl, features=c("MPPED1","FIBCD1","TSPAN18","NECTIN3","AMPH",
                                    "PROX1","SEMA5A",
                                    "GFAP","AQP4","APOE","ID4","SOX2"), 
                   group="domain", swap_rownames = "gene_name",
                   scale=T, center=T, cluster_rows=F, treeheight_col=15, fontsize=7, angle_col=0,
		   fontsize_row=10, fontsize_col=10)

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

sp1 <- plotSpots(sub.manual, annotate="ManualAnnotation", pal=man.pal, point_size=0.4)+
  theme_void()+theme(legend.position="none", aspect.ratio=1)

k17.k18.pal = c("10"="#5ffffb","8"="#99ff99","1"="#61963d", "14"="#add294", "9"="#00dc00",
                "3"="#85FF33", "7"="#00a000", "17"="#B0BF1A", "4"= "#005000", "2"= "#c1c1c1",
                "13"="#444444", "5"="#777777", "15"="#dfa56e",
                "6"="#ff3ffc", "16"="#7a007a", "18"="#ff80fe", "12"="#00006a", "11"="#1e1eff")

sp2 <- plotSpots(sub.manual, annotate="PRECAST_k18", pal=k17.k18.pal, point_size=0.4)+
  theme_void()+theme(legend.position="none", aspect.ratio=1)

ggsave(file="plots/revision/supp_cluster-annotation_neuropil-plots.pdf", 
       gridExtra::grid.arrange(sp1, sp2, hm2[[4]], ncol=3),
       width=7, height=2.5, units="in"
)

############### boxplots for rest of QC metrics
load("plots/spatial_palette_final.rda")
bp1 <- ggplot(as.data.frame(colData(spe)),
             aes(x=domain, y=log2(nuclei_count+1), fill=domain))+
  geom_boxplot(outlier.size = .5, linewidth=.3)+scale_fill_manual(values=spatial.palette)+
  labs(y="VistoSeg nuclei\ncount (log2)")+
  theme_minimal()+theme(legend.position="none", axis.title.x=element_blank(),
                        axis.title.y=element_text(size=8), axis.text=element_text(size=7))

bp2 <- ggplot(as.data.frame(colData(spe)),
              aes(x=domain, y=log10(sum), fill=domain))+
  geom_boxplot(outlier.size = .5, linewidth=.3)+scale_fill_manual(values=spatial.palette)+
  labs(y="Library size\n(log10)")+
  theme_minimal()+theme(legend.position="none", axis.title.x=element_blank(),
                        axis.title.y=element_text(size=8), axis.text=element_text(size=7))

bp3 <- ggplot(as.data.frame(colData(spe)),
              aes(x=domain, y=log10(detected), fill=domain))+
  geom_boxplot(outlier.size = .5, linewidth=.3)+scale_fill_manual(values=spatial.palette)+
  labs(y="Detected genes\n(log10)")+
  theme_minimal()+theme(legend.position="none", axis.title.x=element_blank(),
                        axis.title.y=element_text(size=8), axis.text=element_text(size=7))

bp4 <- ggplot(as.data.frame(colData(spe)),
              aes(x=domain, y=expr_chrM_ratio, fill=domain))+
  geom_boxplot(outlier.size = .5, linewidth=.3)+scale_fill_manual(values=spatial.palette)+
  labs(y="Mitochondrial\nfraction")+
  theme_minimal()+theme(legend.position="none", axis.title.x=element_blank(),
                        axis.title.y=element_text(size=8), axis.text=element_text(size=7))

pdf(file="plots/revision/supp_cluster-annotation_boxplots.pdf", height=5, width=7)
gridExtra::grid.arrange(bp1, bp2, bp3, bp4, ncol=1)
dev.off()
