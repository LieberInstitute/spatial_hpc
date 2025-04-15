# based on: https://github.com/LieberInstitute/spatial_hpc/blob/cebe84643e96a3726b932e3a58c5b4e3ff11f652/code/08_pseudobulk/PRECAST/volcano_plots_by_domain.R
library(dplyr)
library(ggplot2)
library(ggrastr)
library(ggrepel)

load("processed-data/08_pseudobulk/PRECAST/visiumHE_DE_stats_domain.rda")

jtVolcano <- function(domain) {
  fc<-paste0('logFC_',domain)
  fd<-paste0('fdr_',domain)
  tmp = stats$enrichment[,c("gene",fc,fd)]
  colnames(tmp) <- c("gene","fc","fd")
  tmp$sig = tmp$fd<.01 &  abs(tmp$fc)>1
  label.df = filter(tmp, fd<.01, abs(fc)>1)
  lbl1 = slice_max(label.df, n=4, fc)
  lbl2 = slice_min(setdiff(label.df, lbl1), n=8-nrow(lbl1), fd)
  p<- ggplot(tmp, aes(x=fc, y=-log10(fd), color=sig))+
    geom_point(size=.3)+scale_color_manual(values=c("grey","red2"))+
    geom_hline(aes(yintercept=2), lty="longdash", linewidth=.3)+
    geom_vline(aes(xintercept=-1), lty="longdash", linewidth=.3)+geom_vline(aes(xintercept=1), lty="longdash", linewidth=.3)+
    geom_text_repel(data=union(lbl1, lbl2), aes(label=gene), color="black", 
                    min.segment.length = 0, max.overlaps = Inf, size=2.5, fontface="bold")+
    labs(title=domain)+scale_x_continuous(expand=expansion(mult=.2))+
    scale_y_continuous(limits=c(-1,max(-log10(tmp$fd)+5)))+
    theme_minimal()+theme(legend.position="none", axis.line = element_line(color="black", linewidth=.5),
                          #axis.ticks = element_blank(), 
                          panel.grid=element_blank(),
                          axis.title= element_blank(), plot.title.position = "plot",
                          text= element_text(size=10))
  return(rasterize(p,dpi=350,layers = 'Point'))
}

domains<- c('GCL','CA2.4','CA1','SUB','SUB.RHP','RHP','GABA','SL.SR',
            'ML','SR.SLM','SLM.SGZ','WM.1','WM.2','WM.3','Vascular','Choroid') 

pdf(file="plots/revision/volcanoes_for_ED.pdf", height=8, width=8)
do.call(gridExtra::grid.arrange, c(lapply(domains, jtVolcano), ncol=4))
dev.off()
