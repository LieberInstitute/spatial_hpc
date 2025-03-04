library(SpatialExperiment)
library(ggspavis)
library(ggplot2)

load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))
#sub=spe[,spe$brnum=="Br3942"]
sub1 = spe[,spe$sample_id=="V11L05-333_B1"]

load("plots/spatial_palettes.rda")

plotSpots(sub1, annotate="domain", point_size=1)+
  scale_color_manual(values=srt.palette)+ggtitle("title")+
  theme_void()+theme(legend.text=element_text(size=6), plot.title=element_text(hjust=.5, size=18))

plotSpots(sub1, annotate="nmf55", point_size=1)+
  scale_color_gradient(low="grey90", high="black", limits=c(0,max(spe$nmf55)))+
  theme_void()+theme(text=element_text(size=12), plot.title=element_text(hjust=.5, size=18))

sub1$is_gcl = sub1$domain == "GCL"

#crop spot plot by row (height close to top of GCL)
plotSpots(sub1, annotate="array_col", point_size=1)+
  scale_color_viridis_c(option="turbo")

sub2 = sub1[,sub1$array_col<50]

plotVisium(sub2, annotate = "nmf10", highlight = "is_gcl", image = FALSE, spots=T,
           facets=NULL, point_size=.8)+
  scale_color_manual(values=c("grey80","#005000"), guide="none")+guides(color="none")+
  scale_fill_gradient(low="white", high="black", breaks=c(0, 1e-4, 2e-4), 
                      labels=function(x) format(x, scientific=T, digits=1))+
  ggtitle("nmf10")+
  theme(plot.title=element_text(size=18, hjust=.5))

plotVisium(sub2, annotate = "nmf14", highlight = "is_gcl", image = FALSE, spots=T,
           facets=NULL, point_size=.8)+
  scale_color_manual(values=c("grey80","#005000"), guide="none")+guides(color="none")+
  scale_fill_gradient(low="white", high="black", breaks=c(0, 2e-4, 4e-4),
                      labels=function(x) format(x, scientific=T, digits=1))+
  ggtitle("nmf14")+
  theme(plot.title=element_text(size=18, hjust=.5))

ggsave("plots/revision/Figure5_GC-spot-plots.pdf", gridExtra::grid.arrange(p1, p2, ncol=2), 
       bg="white", height=3, width=6, unit="in")
