library(SpatialExperiment)
library(scater)
library(ggspavis)
library(dplyr)
library(ggplot2)
set.seed(123)

load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))

#Figure 2 heatmap
genes2 = c("PPFIA2",
           "KCNG2","KIT","PRKCG",
           "FNDC1","POU3F1","GFRA1",
           "TOX","COL24A1","PART1",
           "TESPA1","KCNH5","MEF2C","GAD2",
           "MIF","APOC1","FABP7","MT1G","MAN1A2","NTRK2","SFRP2",
           "ABCA2","PHLDB1","MTURN","ACTA2","MFAP2",
           "DNAH11")
p1 <- plotGroupedHeatmap(spe, features=genes2, group="domain", swap_rownames="gene_name",
                   center=T, scale=T, zlim=c(-3,4), 
                   assay.type = "logcounts",
                   cluster_rows=F, cluster_cols=F,
                   #colour=viridisLite::magma(n=10))
                   colour=viridisLite::magma(n=30)[c(rep(1,5),2:30)],
                   angle_col=90, silent=T)
ggsave("plots/revision/Figure2_heatmap.pdf", p1, bg="white", width=5, height=5, units="in")
