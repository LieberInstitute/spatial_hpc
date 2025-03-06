library(reactome.db)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

# load NMFs for stats
load(file="snRNAseq_hpc/processed-data/NMF/nmf_final.rda")
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

# prep reactome
xx <- as.list(reactomePATHID2NAME)
reactome.h <- xx[grep("^Homo",xx)]
x <- as.list(reactomePATHID2EXTID)
reactome.h = reactome.h[intersect(names(reactome.h), names(x))]
x.h <- x[names(reactome.h)]
identical(names(x.h), names(reactome.h))
reactome.h = gsub("Homo sapiens: ","",reactome.h)
names(x.h) = reactome.h


#custom plotting function
gseaPlot <- function(nmf, nmf_terms) {
  file1 = grep(nmf, list.files("snRNAseq_hpc/processed-data/revision/"), value=T)
  results = read.csv(paste0("snRNAseq_hpc/processed-data/revision/",file1), row.names=1)
  tmp = arrange(results[results$pathway %in% nmf_terms,], desc(padj))
  #pull rownames of nonzero
  non0.nmf = rownames(loads)[loads[,nmf]>0]
  #convert to entrez id to match with reactome entires
  non0.id = mapIds(org.Hs.eg.db, keys=non0.nmf, keytype="SYMBOL", column="ENTREZID", multiVals = "first")
  #match to gene names
  names(non0.id) = non0.nmf
  #remove gene names without ID
  non0.id = non0.id[!is.na(non0.id)]
  #pull nmf values of mapped genes
  nmf.stats = loads[names(non0.id),nmf]
  #convert names to entrez ID
  names(nmf.stats) = non0.id
  #place in rank order
  nmf.stats = sort(nmf.stats, decreasing=T)
  #make into df
  nmf.df = do.call(rbind, lapply(1:nrow(tmp), function(x) {
    term1 = tmp$pathway[x]
    as.data.frame(list("id"=names(nmf.stats),
                       "weight"=nmf.stats,
                       "rank"=1:length(nmf.stats),
                       "term"=term1,
                       "present"=names(nmf.stats) %in% x.h[[term1]],
                       y=x-1, yend=x))
  }))
  y_labels= sapply(1:nrow(tmp), function(x) paste0("NES= ",round(tmp[x,"NES"],2),
                                                   "\n(",tmp[x,"size"],"/ ",length(x.h[[tmp[x,"pathway"]]]),")\n",
                                                   "padj= ",format(tmp[x,"padj"], scientific=T, digits=2)))
  ggplot(filter(nmf.df, present==T), aes(x=rank, xend=rank, y=y, yend=yend, color=term))+
    geom_segment()+scale_x_continuous(limits=c(0,max(nmf.df$rank)), expand=c(0,0))+
    scale_y_continuous(breaks=c(.5,1.5,2.5), labels=y_labels, expand=c(0,0))+
    theme_minimal()+ggtitle(nmf)+
    theme(axis.title.y=element_blank(), legend.title=element_blank(),
          legend.position="bottom", legend.direction = "vertical",
          panel.grid.minor=element_blank(), panel.grid.major=element_blank())
}




#nmf44.terms = c("L1CAM interactions","Nervous system development","Axon guidance")
#p1 <- gseaPlot("nmf44", nmf44.terms)
nmf77.terms = c("Axon guidance","Signaling by ROBO receptors","Developmental Biology")
p1 <- gseaPlot("nmf77", nmf77.terms)

nmf81.terms = c("Signaling by Receptor Tyrosine Kinases","Extracellular matrix organization","RHO GTPase cycle")
p2 <- gseaPlot("nmf81", nmf81.terms)

nmf13.terms = c("Transmission across Chemical Synapses","Activation of NMDA receptors and postsynaptic events","Protein-protein interactions at synapses")
p3 <- gseaPlot("nmf13", nmf13.terms)

nmf7.terms = c("Transmission across Chemical Synapses","GABA synthesis, release, reuptake and degradation","RHO GTPase Effectors")
p4 <- gseaPlot("nmf7", nmf7.terms)

ggsave("snRNAseq_hpc/plots/revision/Figure4_custom-gsea-plots.pdf", gridExtra::grid.arrange(p1, p2, p3, p4, ncol=1),
       bg="white", height=16, width=12, units="in")
