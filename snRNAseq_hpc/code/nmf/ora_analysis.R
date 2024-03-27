library(ggplot2)
library(dplyr)
library(RcppML)
library(org.Hs.eg.db)
library(clusterProfiler)
library(CoGAPS)

###load sce
load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_nmf_final.rda'))
###load nmf pats
load(file=here::here('snRNAseq_hpc','processed-data','nmf','nmf_final.rda'))

##mark technical patterns for discard
discard<-c('nmf37','nmf28','nmf18','nmf16','nmf94','nmf92','nmf2')

##set up marker gene detection
loads<-x@w
loads<-loads[,!colnames(loads) %in% discard]
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

##now filter mito genes and non-protein coding genes
protein<-rownames(sce)[rowData(sce)$gene_type=='protein_coding']
loads<-loads[rownames(loads) %in% protein,]
mito<-rownames(sce)[which(seqnames(sce) == "chrM")]
loads<-loads[!rownames(loads) %in% mito,]
##now get markers
marks<-patternMarkers(loads,x@h[rownames(x@h) %in% colnames(loads),],'all',1,100)

genes<-marks$PatternMarkers
names(genes)<-rownames(loads)
go<-list()
for(i in 1:length(genes)){
 
go[[i]] <- enrichGO(gene     = chosen,
                            universe      = rownames(loads),
                            OrgDb         = org.Hs.eg.db,
                            #organism='hsa',
                            ont           = "ALL",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.1,
                            readable      = TRUE,
                            keyType = 'SYMBOL')
}

###save go
save(go,file=here::here('snRNAseq_hpc','processed-data','nmf','go_analysis.rda'))