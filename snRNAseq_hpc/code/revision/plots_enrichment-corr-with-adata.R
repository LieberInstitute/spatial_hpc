library(SingleCellExperiment)
library(edgeR)
library(scuttle)
library(spatialLIBD)

reg2 = read.csv("snRNAseq_hpc/python_analysis/processed-data/adata_enrichment_stats_scvi-leiden-k10.csv", row.names=1)
adata_stats<-reg2[,grep("^t_stat",colnames(reg2))]
colnames(adata_stats)<-gsub(pattern='t_stat_X', replacement='c', colnames(adata_stats))
dim(adata_stats)
#21969    34

##############_using superfine labels_################
sn_stats = read.csv("snRNAseq_hpc/processed-data/revision/sn_enrichment_stats_superfine.csv", row.names=1)
sn_stats = sn_stats[,grep("t_stat",colnames(sn_stats))]
colnames(sn_stats)<-gsub(pattern='t_stat_', replacement='', colnames(sn_stats))
dim(sn_stats)
#21104    59

both.genes = intersect(rownames(sn_stats), reg2$gene)
length(both.genes)
#20894
adata_geneid = rownames(reg2)[reg2$gene %in% both.genes]
adata_stats = adata_stats[adata_geneid,]
adata_genes = reg2[rownames(adata_stats),"gene"]
rownames(adata_stats) = adata_genes


sn_stats = sn_stats[both.genes,]

cor.mtx = cor(adata_stats,sn_stats)

c_ordered = c("Astro.2","Astro.3","Astro.1","Oligo.1","Oligo.2","OPC","COP",
              "Micro.1","Micro.2","Macro.Tcell","Ependy",
              "CP.1","CP.2","CP.3","Endo.2","Endo.1","PC.SMC","VLMC",
              "PENK","CORT","SST","PV.FS","C1QL1","CRABP1","LAMP5.MGE","LAMP5.CGE","CXCL14","VIP","HTR3A",
              "Thal",
              "GC.1","GC.2","GC.3","GC.4","GC.5",
              "MC","CA3.1","CA3.2","CA2","CA1","ProS","Sub.1","Sub.2",
              "L2.3.2","L2.3.4","L2.3.3","L2.3.6",
              "L2.3.1","L2.3.5",
              "L5.1","L5.2","L6.1","L6.2","L6b",
              "HATA","AHi.1","AHi.2","AHi.3","AHi.4"
)
r_ordered = c("c28","c4","c2","c15","c32","c8","c10","c12",
              "c25","c5","c31","c17","c18","c19","c23",
              "c0","c26","c20","c3","c27","c1",
              "c24","c11","c30","c21","c13","c29","c14","c33","c6","c16",
              "c9","c7","c22")

pal = c("#440154FF","#440154FF","#440154FF","#440154FF","#482878FF","#3E4A89FF","#31688EFF","#26828EFF","#1F9E89FF","#35B779FF","#6DCD59FF","#B4DE2CFF","#FDE725FF")
pheatmap::pheatmap(cor.mtx[r_ordered,c_ordered],
                   #clustering_method = "ward.D", 
                   cluster_cols=F,cluster_rows=F, 
                   color=pal)
