library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(scater)
library(gridExtra)
set.seed(123)

load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))

#binarize nmfs
spe$sub1 = spe$nmf40>((max(spe$nmf40)*.95)/5)
spe$sub2 = spe$nmf54>((max(spe$nmf54)*.95)/5)
spe$subiculum = spe$sub1|spe$sub2
spe$ca1 = spe$nmf15>((max(spe$nmf15)*.95)/5)
spe$ca3.2 = spe$nmf63>((max(spe$nmf63)*.95)/5)
spe$ca3.1 = spe$nmf11>((max(spe$nmf11)*.95)/5)
spe$ca3 = spe$ca3.1|spe$ca3.2
spe$ENT_sup3 = spe$nmf27>((max(spe$nmf27)*.95)/5)
spe$ENT_sup2 = spe$nmf45>((max(spe$nmf45)*.95)/5)
spe$ENT_sup1 = spe$nmf84>((max(spe$nmf84)*.95)/5)
spe$ENT_sup = spe$ENT_sup1|spe$ENT_sup2|spe$ENT_sup3
spe$ENT_L5 =  spe$nmf51>((max(spe$nmf51)*.95)/5)

spe$binary_label = "other"
colData(spe)[spe$subiculum==T & spe$ca1==F & spe$ca3==F & spe$ENT_sup==F & spe$ENT_L5==F,"binary_label"] = "Subiculum"
colData(spe)[spe$subiculum==F & spe$ca1==T & spe$ca3==F & spe$ENT_sup==F & spe$ENT_L5==F,"binary_label"] = "CA1"
colData(spe)[spe$subiculum==F & spe$ca1==F & spe$ca3==1 & spe$ENT_sup==F & spe$ENT_L5==F,"binary_label"] = "CA3"
colData(spe)[spe$subiculum==F & spe$ca1==F & spe$ca3==F & spe$ENT_sup==T & spe$ENT_L5==F,"binary_label"] = "ENT_sup"
colData(spe)[spe$subiculum==F & spe$ca1==F & spe$ca3==F & spe$ENT_sup==F & spe$ENT_L5==T,"binary_label"] = "ENT_L5"


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



binary.palette = c("CA1"="#984ea3", "CA3"="#fec44f",
                   "Subiculum"="#e41a1c",
                   "ENT_L5"="#377eb8",
                   "ENT_sup"="#61963d",#"L2/3.PrS.Ent"="#add294",
                   "other"="grey85")

brains <- c("Br3942","Br6522","Br8667","Br2743","Br6423","Br2720","Br6432","Br8492","Br6471","Br8325")
p<-list()
for (j in seq_along(brains)) {
  speb <- spe[, (colData(spe)$brnum == brains[j])]
  speb$sample_id <- droplevels(speb$sample_id)
  p[[j]]<-plotVisium(speb, spots=T, image=F, highlight="binary_label", annotate="nmf17",
                     facets="sample_id",point_size=.5)+ 
    facet_grid(rows=vars(facet_row), cols=vars(facet_column))+
    ggtitle(brains[j])+scale_y_continuous(expand=expansion(add=4))+scale_x_continuous(expand=expansion(add=4))+
    scale_color_manual(values=binary.palette, guide="none")+
    scale_fill_gradient(low="white", high="black", limits=c(0, max(spe$nmf17)), guide="none")+
    guides(color=NULL, fill=NULL)+
    #theme_bw()+
    theme(#legend.text = element_text(size = 8),
      plot.title = element_text(size = 16, hjust=.5, margin=unit(c(5,0,0,0), "pt")),
      strip.text=element_blank(),
      panel.spacing=unit(5,'points'),plot.margin=unit(c(0,5,0,0), "pt"),
    )
}

lplot = ggplot(cbind.data.frame("y"=seq(binary.palette), "labels"=names(binary.palette),
                                 "nmf17"=seq(from=0, to=max(spe$nmf17), length.out=length(binary.palette))),
               aes(x=1, y=y, color=labels, fill=nmf17))+
  geom_point(size=3)+scale_color_manual(values=binary.palette, guide="none")+
  geom_text(aes(label=labels), nudge_x = .02, hjust=0, color="black", size=4)+
    scale_fill_gradient(low="white", high="black", 
                        labels=function(x) format(x, scientific=T, digits=1),
                        guide=guide_colorbar(theme=theme(legend.position="bottom", 
                                                              legend.direction="vertical")))+
  scale_y_reverse()+scale_x_continuous(expand=expansion(add=c(.1,.2)))+#theme_bw()+
  theme_void()+
  theme(legend.position="bottom", #legend.justification="center",
                     plot.margin=unit(c(2,.5,2,.5),"cm"))

laymat=rbind(c(1,1,2,2,3,3),c(1,1,2,2,3,3),
             c(4,4,5,5,6,6),c(4,4,5,5,6,6),
             c(7,NA,8,10,10,11),c(7,NA,8,10,10,11),
             c(9,9,9,10,10,11))

pdf(file = "plots/revision/all-brain_spotplot_binary-nmf.pdf",
    width=9, height=12)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],lplot, 
             layout_matrix=laymat)
dev.off()
