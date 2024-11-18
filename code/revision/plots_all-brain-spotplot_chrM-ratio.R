library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(scater)
library(gridExtra)
set.seed(123)

load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))

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
  p[[j]]<-plotVisium(speb, spots=T, image=F, annotate="expr_chrM_ratio", #highlight="expr_chrM_ratio",
                      facets="sample_id",point_size=.6)+ 
    facet_grid(rows=vars(facet_row), cols=vars(facet_column))+
    ggtitle(brains[j])+scale_y_continuous(expand=expansion(add=4))+scale_x_continuous(expand=expansion(add=4))+
    scale_fill_gradient(low="grey90",high="black", limits=c(0,.6), guide="none")+
    #scale_color_manual(values="transparent", aesthetics="highlight")+
    guides(color=NULL, fill=NULL)+
    #theme_bw()+
    theme(#legend.text = element_text(size = 8),
      plot.title = element_text(size = 16, hjust=.5, margin=unit(c(5,0,0,0), "pt")),
      strip.text=element_blank(),
      panel.spacing=unit(5,'points'),plot.margin=unit(c(0,5,0,0), "pt"),
    )
}

lplot = ggplot(cbind.data.frame("y"=c(1:10), "labels"=c(1,2),
                                "chrM.ratio"=seq(from=0,to=.6, length.out=10)),
               aes(x=1, y=y, fill=chrM.ratio))+
  geom_tile(stat="identity", alpha=0)+
  scale_fill_gradient(low="grey90", high="black", 
                      guide=guide_colorbar())+theme_void()+theme(legend.position="inside")


laymat=rbind(c(1,1,2,2,3,3),c(1,1,2,2,3,3),
             c(4,4,5,5,6,6),c(4,4,5,5,6,6),
             c(7,NA,8,10,10,11),c(7,NA,8,10,10,11),
             c(9,9,9,10,10,11))

pdf(file = "plots/revision/all-brain_spotplot_chrM-ratio.pdf",
    width=9, height=12)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],lplot, 
             layout_matrix=laymat)
dev.off()
