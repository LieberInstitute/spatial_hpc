library(SpatialExperiment)
library(SingleCellExperiment)
library(scran)
library(dplyr)
library(ggplot2)
library(scater)
set.seed(123)

load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))

#nrn heatmap - broad strokes markers
glia_clus = c(34,36,62,#CP
              19,47,41,30,#endo
              48,#ependy
              8,20,27,#astro
              1,9,#oligo
              44,24,#opc
              2,6,42)#microglia
nrn_clus = setdiff(unique(sce$k_5_louvain), glia_clus)
sce_nrn = sce[,sce$k_5_louvain %in% nrn_clus]

partial.exc = c(59,18,39,40,29,#4,
                3,52,56,10,16,7,
                58,#cajal
                33,55,49,12,31,
                14,26,
                45,11,38,
                21,61,
                25,
                22,32,35,15,
                5,43,
                13,51,28,
                54,37,
                50,
                4,#thal
                #61,
                57,60,53,46#HAT
)

sce_nrn$k_5_louvain = factor(sce_nrn$k_5_louvain, levels=partial.exc) 
plotGroupedHeatmap(sce_nrn, features=c("GAD2","ADARB2","LHX6",#"LAMP5",
                                       "RELN",#"PENK","MEIS2",
                                       #"SLC17A7",
                                       "PROX1","CALB1",
                                       "FIBCD1","TSPAN18","CARTPT","FN1",
                                       "TLE4","SATB2","CUX2"), 
                   group="k_5_louvain", cluster_rows=F, cluster_cols=F,
                   scale=T, center=T, zlim=c(-2,4), angle_col=0)#+

#thalamus boxplot
thal.plot = mutate(as.data.frame(colData(sce_nrn)), tcf7l2=logcounts(sce_nrn)["TCF7L2",],
                   is_thal=k_5_louvain==4)
ggplot(thal.plot, aes(x=k_5_louvain, y=tcf7l2, fill=is_thal))+
  geom_boxplot(outlier.size=.3)+scale_fill_manual(values=c("grey","skyblue"))+
  theme_bw()


#annotate amygdala by looking at SPE DEGs for amygdala 
load(here::here('processed-data','06_clustering','PRECAST','spe_precast_HE_domain.rda'))
table(spe$VSPG)

#from https://github.com/LieberInstitute/spatial_hpc/blob/main/code/08_pseudobulk/PRECAST/create_pseudobulk_spe_captureArea.R
spe$amy<-ifelse(spe$sample_id %in%
                         c('Br6423_V10B01-085_C1','Br6432_V10B01-086_B1') &
                         spe$domain %in% c('RHP'),T,F)
spe$amy<-ifelse(spe$sample_id =='Br6423_V10B01-085_A1',T,spe$amy) #one sample is entirely amygdala
#6423, 6432, 8667

table(spe$amy)
#FALSE   TRUE 
#146306   4611 

spe_sce<-SingleCellExperiment(assays=list(counts=counts(spe), logcounts=logcounts(spe)),
                              colData=colData(spe),rowData=rowData(spe))

#amy <- spe_sce[,spe_sce$brnum %in% amy.brnum & spe_sce$domain=="RHP" & spe_sce$sample_id!='Br6423_V10B01-085_A1']
amy <- spe_sce[,spe_sce$domain=="RHP"]
table(amy$amy)
#FALSE  TRUE 
#3490  4519 
amy.markers <-findMarkers(amy, groups=amy$amy, pval.type='all', direction='up')
amy.markers.up = mutate(as.data.frame(amy.markers[[2]]), gene_name=rowData(amy)[rownames(amy.markers[[2]]),"gene_name"]) %>%
  filter(FDR<.05 & summary.logFC>.3)
nrow(amy.markers.up)

amy.markers.dn = mutate(as.data.frame(amy.markers[[1]]), gene_name=rowData(amy)[rownames(amy.markers[[1]]),"gene_name"]) %>%
  filter(FDR<.05 & summary.logFC>.3)
nrow(amy.markers.dn)

#now that we have markers, subset to just clusters with strong cortical layer markers
cort = c("Sub.2","L6.1","L6.2","L6b","L5.2","L5.1","L2/3.3","L2/3.2","L2/3.1","L2/3.6","L2/3.4","L2/3.5","HATA","AHi.1","AHi.2","AHi.3","AHi.4")
sce_cort = sce[,sce$superfine.cell.class %in% cort]
sce_cort$superfine.cell.class<-droplevels(sce_cort$superfine.cell.class)
dim(sce_cort)

#plot up and down amy markers to get row order for final plot
p1 = plotGroupedHeatmap(sce_cort, features=amy.markers.up$gene_name,
                   group="superfine.cell.class",
                   scale=T, center=T)
up_order = p1$tree_row$labels[p1$tree_row$order]
p2 <- plotGroupedHeatmap(sce_cort, features=amy.markers.dn$gene_name, 
                   group="superfine.cell.class",
                   scale=T, center=T)
dn_order = p2$tree_row$labels[p2$tree_row$order]
clus_order = p2$tree_col$labels[p2$tree_col$order]
sce_cort$superfine.cell.class = factor(sce_cort$superfine.cell.class, levels=clus_order)
sce_cort$k_5_louvain = factor(sce_cort$k_5_louvain, levels=c(28,50,54,13,
                                                             15,22,5,
                                                             61,46,57,53,60, #amyg
                                                             32,35,51,43,37))
plotGroupedHeatmap(sce_cort, features=c(up_order, dn_order),
                   group="k_5_louvain", border_color=NA,
                   scale=T, center=T, cluster_rows=F, cluster_cols=F
)
