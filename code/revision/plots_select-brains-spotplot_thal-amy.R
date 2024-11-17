library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(scater)
library(gridExtra)
set.seed(123)

load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))
load("plots/spatial_palette_final.rda")

##################### thalamus
spe.thal = spe[,spe$brnum=="Br8325"]
spe.thal$sample_id = factor(spe.thal$sample_id, levels=c("V11A20-297_C1","V11A20-297_D1","V11A20-297_A1","V11L05-335_D1","V11A20-297_B1"))
spe.thal$facet_column = factor(spe.thal$sample_id, levels=levels(spe.thal$sample_id),
                               labels=c("c1","c2","c1","c2","c2"))
spe.thal$facet_row = factor(spe.thal$sample_id, levels=levels(spe.thal$sample_id),
                               labels=c("r1","r1","r2","r2","r3"))
  
pt1 <- plotVisium(spe.thal, spots=T, image=F, highlight="domain", annotate="SHOX2",
                  facets="sample_id", point_size=.8, assay="logcounts")+ 
  facet_grid(rows=vars(facet_row), cols=vars(facet_column))+
  scale_y_continuous(expand=expansion(add=4))+scale_x_continuous(expand=expansion(add=4))+
  scale_color_manual(values=spatial.palette, guide="none")+
  scale_fill_gradient(low="white", high="black", guide="none")+#, limits=c(0, max(assay(spe, "counts1")["TCF7L2",])), guide="none")+
  guides(color=NULL, fill=NULL)+
  theme(strip.text=element_blank(),
        panel.spacing=unit(5,'points'),plot.margin=unit(c(0,0,5,0), "pt"))

lplot = ggplot(cbind.data.frame("y"=seq(spatial.palette), "labels"=names(spatial.palette),
                                "SHOX2"=seq(from=0, to=max(logcounts(spe.thal)["SHOX2",]), 
                                              length.out=length(spatial.palette))),
               aes(x=1, y=y, color=labels, fill=SHOX2))+
  geom_point(size=3)+scale_color_manual(values=spatial.palette, guide="none")+
  geom_text(aes(label=labels), nudge_x = .06, hjust=0, color="black", size=4)+
  scale_fill_gradient(low="white", high="black", 
                      #labels=function(x) format(x, scientific=T, digits=1),
                      guide=guide_colorbar(theme=theme(legend.position="right", 
                                                       legend.direction="vertical")))+
  guides(color=NULL)+
  scale_y_reverse()+scale_x_continuous(expand=expansion(add=c(.1,.4)))+
  theme_void()+theme(plot.margin=unit(c(2,0,2,.5),"cm"))

laymat = rbind(c(1,1,2),c(1,1,2),c(1,1,2))
grid.arrange(pt1,lplot, layout_matrix=laymat, top="Br8325")

#H&E image
plotVisium(spe.thal, spots=F, image=T, #highlight="domain", annotate="SHOX2",
           facets="sample_id", point_size=.8)#+ 
  #facet_grid(rows=vars(facet_row), cols=vars(facet_column))+
  #scale_y_continuous(expand=expansion(add=4))+scale_x_continuous(expand=expansion(add=4))+
  #scale_color_manual(values=spatial.palette, guide="none")+
  #scale_fill_gradient(low="white", high="black", guide="none")+#, limits=c(0, max(assay(spe, "counts1")["TCF7L2",])), guide="none")+
  #guides(color=NULL, fill=NULL)+
  theme(strip.text=element_blank(),
        panel.spacing=unit(5,'points'),plot.margin=unit(c(0,0,5,0), "pt"))

##################### amygdala
spe.amy = spe[,spe$brnum=="Br6423"]
spe.amy$sample_id = factor(spe.amy$sample_id, levels=c("V10B01-085_B1","V10B01-085_A1","V10B01-085_D1","V10B01-085_C1"))

pa<-list()
for (j in seq_along(levels(spe.amy$sample_id))) {
  spea <- spe.amy[,spe.amy$sample_id == levels(spe.amy$sample_id)[[j]]]
  pa[[j]]<-plotVisium(spea, spots=T, image=F, highlight="domain", annotate="SLC17A6",
                     point_size=.8, assay="logcounts")+ 
    scale_y_continuous(expand=expansion(add=4))+scale_x_continuous(expand=expansion(add=4))+
    scale_color_manual(values=spatial.palette, guide="none")+
    scale_fill_gradient(low="white", high="black", limits=c(0, max(logcounts(spe.amy)["SLC17A6",])), 
                        guide="none")+
    guides(color=NULL, fill=NULL)+theme(strip.text=element_blank())
}  

lplot = ggplot(cbind.data.frame("y"=seq(spatial.palette), "labels"=names(spatial.palette),
                                "SLC17A6"=seq(from=0, to=max(logcounts(spe.amy)["SLC17A6",]), 
                                             length.out=length(spatial.palette))),
               aes(x=1, y=y, color=labels, fill=SLC17A6))+
  geom_point(size=3)+scale_color_manual(values=spatial.palette, guide="none")+
  geom_text(aes(label=labels), nudge_x = .06, hjust=0, color="black", size=4)+
  scale_fill_gradient(low="white", high="black", 
                      #labels=function(x) format(x, scientific=T, digits=1),
                      guide=guide_colorbar(theme=theme(legend.position="right", 
                                                       legend.direction="vertical")))+
  guides(color=NULL)+
  scale_y_reverse()+scale_x_continuous(expand=expansion(add=c(.1,.4)))+
  theme_void()+theme(plot.margin=unit(c(2,0,2,.5),"cm"))

laymat = rbind(c(NA,2,5),c(1,2,5),c(1,2,5),c(1,2,5),c(1,2,5),
               c(3,4,5),c(3,4,5),c(3,4,5),c(3,4,5))

grid.arrange(pa[[1]],pa[[2]],pa[[3]],pa[[4]],lplot, layout_matrix=laymat, top="Br6423")

#H&E images
plotVisium(spe.amy, spots=F, image=T, #highlight="domain", annotate="SHOX2",
           facets="sample_id", point_size=.8)#+ 
