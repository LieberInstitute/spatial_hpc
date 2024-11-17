library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(scater)
library(gridExtra)
set.seed(123)

load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))


spe$is_rhp = spe$domain=="RHP"

#below is necessary because unlike plotSpots, there is a bug in plotVisium where the sanity check
#for the assay name of 'length(grep(assay, assayNames(spe))) == 1' must be true but because 
#"counts" and "logcounts" both have the grep, the function won't proceed.
#I tried setting the assay to "^counts", and this bypasses the sanity check but the actual 
#fetch assay call requires the true assay name so this errors as well:
### 'assay(<SpatialExperiment>, i="character", ...)' invalid subscript 'i'
### '^counts' not in names(assays(<SpatialExperiment>))
assayNames(spe) <- c("counts1","logcounts")

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

#broad.palette = c("Neuron"="#00a000", "Neuropil"="#dfa56e", "WM"="#ff80fe", "Vasc_CSF"="#00006a")
broad.palette2= c("Neuron"="#add294", "Neuropil"="#efd2b6", "WM"="#ffccfe", "Vasc_CSF"="#9999c3")

brains <- c("Br3942","Br6522","Br8667","Br2743","Br6423","Br2720","Br6432","Br8492","Br6471","Br8325")
p<-list()
for (j in seq_along(brains)) {
  speb <- spe[, (colData(spe)$brnum == brains[j])]
  speb$sample_id <- droplevels(speb$sample_id)
  p[[j]]<-plotVisium(speb, spots=T, image=F, highlight="broad.domain", annotate="TCF7L2",
                     facets="sample_id",point_size=.5, assay="counts1")+ 
    facet_grid(rows=vars(facet_row), cols=vars(facet_column))+
    ggtitle(brains[j])+scale_y_continuous(expand=expansion(add=4))+
    scale_x_continuous(expand=expansion(add=4))+
    scale_color_manual(values=broad.palette2, guide="none")+
    scale_fill_gradient(low="white", high="black", limits=c(0, max(assay(spe, "counts1")["TCF7L2",])), guide="none")+
    guides(color=NULL, fill=NULL)+
    #theme_bw()+
    theme(#legend.text = element_text(size = 8),
      plot.title = element_text(size = 16, hjust=.5, margin=unit(c(5,0,0,0), "pt")),
      strip.text=element_blank(),
      panel.spacing=unit(5,'points'),plot.margin=unit(c(0,5,0,0), "pt"),
    )
}

lplot = ggplot(cbind.data.frame("y"=seq(broad.palette2), "labels"=names(broad.palette2),
                                "TCF7L2"=seq(from=0, to=max(assay(spe, "counts1")["TCF7L2",]), 
                                            length.out=length(broad.palette2))),
               aes(x=1, y=y, color=labels, fill=TCF7L2))+
  geom_point(size=3)+scale_color_manual(values=broad.palette2, guide="none")+
  geom_text(aes(label=labels), nudge_x = .02, hjust=0, color="black", size=4)+
  scale_fill_gradient(low="white", high="black", 
                      #labels=function(x) format(x, scientific=T, digits=1),
                      guide=guide_colorbar(theme=theme(legend.position="bottom", 
                                                       legend.direction="vertical")))+
  scale_y_reverse(expand=expansion(add=c(.5,.5)))+scale_x_continuous(expand=expansion(add=c(.1,.2)))+#theme_bw()+
  theme_void()+
  theme(legend.position="bottom", #legend.justification="center",
        plot.margin=unit(c(2,.5,2,.5),"cm"))

laymat=rbind(c(1,1,2,2,3,3),c(1,1,2,2,3,3),
             c(4,4,5,5,6,6),c(4,4,5,5,6,6),
             c(7,NA,8,10,10,11),c(7,NA,8,10,10,11),
             c(9,9,9,10,10,11))

pdf(file = "plots/revision/all-brain_spotplot_broad-domain-TCF7L2.pdf",
    width=9, height=12)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],lplot, 
             layout_matrix=laymat)
dev.off()





brains <- c("Br3942","Br6522","Br8667","Br2743","Br6423","Br2720","Br6432","Br8492","Br6471","Br8325")
p<-list()
for (j in seq_along(brains)) {
  speb <- spe[, (colData(spe)$brnum == brains[j])]
  speb$sample_id <- droplevels(speb$sample_id)
  p[[j]]<-plotVisium(speb, spots=T, image=F, highlight="broad.domain", annotate="OPRM1",
                     facets="sample_id",point_size=.5, assay="counts1")+ 
    facet_grid(rows=vars(facet_row), cols=vars(facet_column))+
    ggtitle(brains[j])+scale_y_continuous(expand=expansion(add=4))+
    scale_x_continuous(expand=expansion(add=4))+
    scale_color_manual(values=broad.palette2, guide="none")+
    scale_fill_gradient(low="white", high="black", limits=c(0, max(assay(spe, "counts1")["OPRM1",])), guide="none")+
    guides(color=NULL, fill=NULL)+
    #theme_bw()+
    theme(#legend.text = element_text(size = 8),
      plot.title = element_text(size = 16, hjust=.5, margin=unit(c(5,0,0,0), "pt")),
      strip.text=element_blank(),
      panel.spacing=unit(5,'points'),plot.margin=unit(c(0,5,0,0), "pt"),
    )
}

lplot = ggplot(cbind.data.frame("y"=seq(broad.palette2), "labels"=names(broad.palette2),
                                "OPRM1"=seq(from=0, to=max(assay(spe, "counts1")["OPRM1",]), 
                                            length.out=length(broad.palette2))),
               aes(x=1, y=y, color=labels, fill=OPRM1))+
  geom_point(size=3)+scale_color_manual(values=broad.palette2, guide="none")+
  geom_text(aes(label=labels), nudge_x = .02, hjust=0, color="black", size=4)+
  scale_fill_gradient(low="white", high="black", 
                      #labels=function(x) format(x, scientific=T, digits=1),
                      guide=guide_colorbar(theme=theme(legend.position="bottom", 
                                                       legend.direction="vertical")))+
  scale_y_reverse(expand=expansion(add=c(.5,.5)))+scale_x_continuous(expand=expansion(add=c(.1,.2)))+#theme_bw()+
  theme_void()+
  theme(legend.position="bottom", #legend.justification="center",
        plot.margin=unit(c(2,.5,2,.5),"cm"))

laymat=rbind(c(1,1,2,2,3,3),c(1,1,2,2,3,3),
             c(4,4,5,5,6,6),c(4,4,5,5,6,6),
             c(7,NA,8,10,10,11),c(7,NA,8,10,10,11),
             c(9,9,9,10,10,11))

pdf(file = "plots/revision/all-brain_spotplot_broad-domain-OPRM1.pdf",
    width=9, height=12)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],lplot, 
             layout_matrix=laymat)
dev.off()
