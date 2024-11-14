#library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(pheatmap)
#library(scater)
#library(edgeR)
#library(scuttle)
#library(spatialLIBD)
set.seed(123)

reg <- read.csv("processed-data/revision/spe_enrichment_stats_cluster-k18-plus-amy-thal.csv", row.names=1)
spe_stats = reg[,grep("t_stat",colnames(reg))]
colnames(spe_stats)<-gsub(pattern='t_stat_',
                          replacement='',
                          colnames(spe_stats))
spe_stats = reg[,grep("logFC_",colnames(reg))]
colnames(spe_stats)<-gsub(pattern='logFC_',
                          replacement='',
                          colnames(spe_stats))

sn_stats = read.csv("snRNAseq_hpc/processed-data/revision/sn_enrichment_stats_superfine.csv", row.names=1)
sn_stats = sn_stats[,grep("t_stat",colnames(sn_stats))]
colnames(sn_stats)<-gsub(pattern='t_stat_',
                            replacement='',
                            colnames(sn_stats))
sn_stats = sn_stats[,grep("logFC",colnames(sn_stats))]
colnames(sn_stats)<-gsub(pattern='logFC_',
                         replacement='',
                         colnames(sn_stats))

replace.names = reg$gene#rowData(spe_pseudo)[rownames(spe_stats),"gene_name"]
grep("HSPA14",replace.names)
#domain
replace.names[6527] = "HSPA14_ENSG00000284024"
replace.names[6528] = "HSPA14_ENSG00000187522"
#cluster
replace.names[6556] = "HSPA14_ENSG00000284024"
replace.names[6557] = "HSPA14_ENSG00000187522"

rownames(spe_stats) = replace.names

new.names = intersect(replace.names, rownames(sn_stats))
length(new.names)
#13364 -- domain
#13395 -- cluster

spe_stats = spe_stats[new.names,]
sn_stats = sn_stats[new.names,]
dim(sn_stats)

cor.mtx = cor(spe_stats,sn_stats)

pal = c("#440154FF","#440154FF","#440154FF","#440154FF","#482878FF","#3E4A89FF","#31688EFF","#26828EFF","#1F9E89FF","#35B779FF","#6DCD59FF","#B4DE2CFF","#FDE725FF")
pheatmap(cor.mtx[r_ordered,c_ordered],
                   #clustering_method = "ward.D", 
                   cluster_cols=F,cluster_rows=F, 
                   color=pal)
#                   color=viridis::viridis(100))

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
r_ordered = c("SLM.SGZ","SR.SLM","WM.1","WM.2","WM.3","Choroid","Vascular",
              "SL.SR",
              "GABA","Thalamus",
              "ML","GCL",
              "CA2.4.2","CA2.4.1","CA1.1","CA1.2",
              "SUB","SUB.RHP","RHP","Amygdala"
              )
