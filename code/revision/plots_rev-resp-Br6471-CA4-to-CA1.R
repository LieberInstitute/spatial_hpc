library(SingleCellExperiment)
library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(scater)
library(ggspavis)

load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))
load("plots/spatial_palette_final.rda")

#spe$mc.nmf52 = spe$nmf52>((max(spe$nmf52)*.95)/5)
#spe$ca2.nmf61 = spe$nmf61>((max(spe$nmf61)*.95)/5)

#Br3942, all domains, CA1, ProS, Sub1, Sub2
sub1 = spe[,spe$brnum=="Br6471"] #
sub1$sample_id = factor(sub1$sample_id, levels=c("V11L05-335_A1","V11L05-335_C1","V11L05-335_B1"))
sub1$facet_column = factor(sub1$sample_id, levels=levels(sub1$sample_id), labels=c("c1","c2","c3"))

p1 <- plotVisium(sub1, spots=T, annotate="nmf52", highlight="domain",
                 facets="sample_id", image = FALSE, point_size=.8)+coord_cartesian()+
  facet_grid(cols=vars(facet_column), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf52)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  ggtitle("MC/CA4 (nmf52)")+
  theme(legend.text = element_text(size = 8), 
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(3,5,3,0), "pt"))

p2 <- plotVisium(sub1, spots=T, annotate="nmf11", highlight="domain",
                 facets="sample_id", image = FALSE, point_size=.8)+coord_cartesian()+
  facet_grid(cols=vars(facet_column), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf11)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  ggtitle("CA3 (nmf11)")+
  theme(legend.text = element_text(size = 8), 
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(3,5,3,0), "pt"))

p3 <- plotVisium(sub1, spots=T, annotate="nmf63", highlight="domain",
                 facets="sample_id", image = FALSE, point_size=.8)+coord_cartesian()+
  facet_grid(cols=vars(facet_column), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf63)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  ggtitle("CA3 (nmf63)")+
  theme(legend.text = element_text(size = 8), 
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(3,5,3,0), "pt"))

p4 <- plotVisium(sub1, spots=T, annotate="nmf61", highlight="domain",
                 facets="sample_id", image = FALSE, point_size=.8)+coord_cartesian()+
  facet_grid(cols=vars(facet_column), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf61)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  ggtitle("CA2 (nmf61)")+
  theme(legend.text = element_text(size = 8), 
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(3,5,3,0), "pt"))

p5 <- plotVisium(sub1, spots=T, annotate="nmf15", highlight="domain",
                 facets="sample_id", image = FALSE, point_size=.8)+coord_cartesian()+
  facet_grid(cols=vars(facet_column), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf15)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  ggtitle("CA1 (nmf15)")+
  theme(legend.text = element_text(size = 8), 
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(3,5,3,0), "pt"))

pdf(file="plots/revision/rev-resp_Br6471_CA4-to-CA1.pdf", width=7, height=11)
gridExtra::grid.arrange(p1, p2, p3, p4, p5, ncol=1)
dev.off()
