library(SingleCellExperiment)
library(dplyr)
library(ggplot2)

set.seed(123)
#load sce
load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))
avg.expr = rowMeans(logcounts(sce))
#load nmfs
load(file=here::here('snRNAseq_hpc','processed-data','NMF','nmf_final.rda'))
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]
avg.expr = avg.expr[rownames(loads)]

df1 = cbind.data.frame(loads, avg.expr)

plist1 <- lapply(colnames(loads), function(x) {
  tmp = df1[,c(x,"avg.expr")]
  colnames(tmp) = c("nmf","avg.expr")
  tmp$rank = (1+nrow(loads))-rank(tmp$nmf)
  tmp$top150 = tmp$rank<=10
  ggplot(tmp, aes(x=avg.expr, y=nmf, color=top150))+
    geom_point(size=.3)+theme_bw()+
    geom_point(data=filter(tmp, top150==T), size=.5)+
    scale_color_manual(values=c("black","red"))+
    labs(title=x, y="gene weights", x="avg. expr")+
    theme(text=element_text(size=8), #axis.title.y=element_blank()), 
          axis.text=element_blank(), legend.position="none")
})

ggsave(filename="snRNAseq_hpc/plots/revision/nmf-weight_gene-expr_scatter_top10-red-enlarge.png", 
       plot=do.call(gridExtra::grid.arrange, c(plist1, ncol=10)), 
       width=14, height=14, units="in", bg="white")
