library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scater)
library(org.Hs.eg.db)
set.seed(123)

#load sce
load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))

#load nmfs
load(file=here::here('snRNAseq_hpc','processed-data','NMF','nmf_final.rda'))
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

#check for sex
xist = sapply(colnames(loads), function(x) {
  tmp = names(sort(loads[,x], decreasing=T)[1:50])
  "XIST" %in% tmp
})
xist[xist] #just 1, nmf37


chr.loc <- org.Hs.egCHR
ezid = mapIds(org.Hs.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENTREZID", multiVals = "first")
names(ezid) = rownames(sce)
ezid2 = ezid[!is.na(ezid)]
test = as.list(chr.loc[ezid2])
chrY = ezid2[grep('Y',test)]

chrYgenes = sapply(colnames(loads), function(x) {
  tmp = names(sort(loads[,x], decreasing=T)[1:50])
  names(chrY) %in% tmp
})
colSums(chrYgenes) #just 1 pattern has more than 1 in top 50: nmf28

p1 <- plotReducedDim(sce, dimred="UMAP", color_by="nmf37", point_size=.1)+
  scale_color_viridis_c(option="F", direction=-1)+
  labs(title="nmf37", color="nuclei\nweight")+theme_bw()+
  theme(legend.box.spacing = unit(0,"pt"), legend.box.margin = margin(0,0,0,0,"pt"),
        axis.title=element_blank(), legend.key.size=unit(10,"pt"), legend.text=element_text(size=8, margin=margin(0,0,0,3,"pt")))
p2 <- plotReducedDim(sce, dimred="UMAP", color_by="XIST", point_size=.1, by_exprs_values = "logcounts")+
  scale_color_gradient(low="grey90", high="black")+
  labs(title="XIST", color="log2\nCPM")+theme_bw()+
  theme(legend.box.spacing = unit(0,"pt"), legend.box.margin = margin(0,0,0,0,"pt"),
        axis.title=element_blank(), legend.key.size=unit(10,"pt"))
p3 <- plotReducedDim(sce, dimred="UMAP", color_by="nmf28", point_size=.1)+
  scale_color_viridis_c(option="F", direction=-1)+
  labs(title="nmf28", color="nuclei\nweight")+theme_bw()+
  theme(legend.box.spacing = unit(0,"pt"), legend.box.margin = margin(0,0,0,0,"pt"),
        axis.title=element_blank(), legend.key.size=unit(10,"pt"), legend.text=element_text(size=8, margin=margin(0,0,0,3,"pt")))
p4 <- plotReducedDim(sce, dimred="UMAP", color_by="USP9Y", point_size=.1, by_exprs_values = "logcounts")+
  scale_color_gradient(low="grey90", high="black")+
  labs(title="USP9Y", color="log2\nCPM")+theme_bw()+
  theme(legend.box.spacing = unit(0,"pt"), legend.box.margin = margin(0,0,0,0,"pt"),
        axis.title=element_blank(), legend.key.size=unit(10,"pt"))

ggsave(file="snRNAseq_hpc/plots/revision/supp_nmf-sex-specific.png", bg="white",
       gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2), height=9, width=8)
