####project human NMF patterns to mouse snRNAseq ECS data
library(SingleCellExperiment)
library(scater)
library(scran)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RcppML)
library(pheatmap)
library(ggrastr)

####load NMF patterns
load(file=here::here('snRNAseq_hpc','processed-data','NMF','nmf_final.rda'))

###load mouse sce (sce.subset)
load(file=here::here('snRNAseq_hpc','processed-data','NMF','sce_subset.rda'))


# Translate from one species to the other using the orthology
orthology<-read.csv(file=here::here('snRNAseq_hpc','processed-data',
                                    'NMF','human_mouse_orthologs.csv'))
names <- orthology[orthology$Column3 %in% rownames(sce.subset),]

names <- names[match(rownames(sce.subset), names$Column3),]

setdiff(names$Column3, rownames(sce.subset))

rownames(sce.subset) <- names$Column1

########projection#####
# Project human snRNAseq NMF patterns onto mouse snRNAseq

set.seed(1029)
i<-intersect(rownames(sce.subset),rownames(x@w))
loadings<-x@w
loadings<-loadings[rownames(loadings) %in% i,]
sce.subset2<-sce.subset[rownames(sce.subset) %in% i,]
loadings<-loadings[match(rownames(sce.subset2),rownames(loadings)),]
proj<-project(loadings,logcounts(sce.subset2),L1=0)

##rescale patterns and add to sce
proj<-t(proj)
proj<-apply(proj,2,function(x){x/sum(x)})

###add to sce
colData(sce.subset)<-cbind(colData(sce.subset),proj)

#########supp figure plots###########
###correlation with seizures
##first rename the levels of condition
levels(sce.subset$condition)<-c('control','activated')

##now make one-hot encoded matrix
library(reshape2)
data<-as.data.frame(sce.subset$condition)
colnames(data)<-'condition'
onehot_condition <-  dcast(data = data, rownames(data) ~ condition, length)
rownames(onehot_condition)<-onehot_condition[,1]
onehot_condition[,1]<-as.numeric(onehot_condition[,1])
onehot_condition<-onehot_condition[order(onehot_condition[,1],decreasing=F),]
onehot_condition[,1]<-NULL

##correlate
cors<-cor(onehot_condition,proj)
cors<-cors[,-2] ##nmf2 is not expressed at all, comes up as NA. need to remove

##order cors by "activated
cors<-cors[,order(cors[2,],decreasing=T)]

##make heatmap
pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','heatmap_ecs_nmf.pdf'),h=4,w=8.5)
pheatmap(cors,cluster_rows=F,cluster_cols=F)
dev.off()


####violin plots
sce.pats<-SingleCellExperiment(assays=list(counts=t(proj)))
logcounts(sce.pats)<-counts(sce.pats)
colData(sce.pats)<-colData(sce.subset)

rownames(sce.pats)[match(c(
                        'nmf5','nmf26','nmf10','nmf14',
                        'nmf52','nmf11','nmf63','nmf61','nmf15',
                       'nmf32','nmf40','nmf54','nmf22','nmf65',
                       'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                       'nmf84','nmf62','nmf69',
                       'nmf83','nmf73','nmf50','nmf67','nmf35','nmf47',
                       'nmf88','nmf46','nmf74','nmf60'
                       ),rownames(sce.pats))]<-c(
                           'nmf5 GC','nmf26 GC','nmf10 GC','nmf14 GC',
                           'nmf52 MC','nmf11 CA3.1','nmf63 CA3.2','nmf61 CA2','nmf15 CA1',
                           'nmf32 ProS','nmf40 Sub.1','nmf54 Sub.2','nmf22 L6.2','nmf65 L6/6b',
                           'nmf53 L6b','nmf51 L5/6.1','nmf68 L5/6','nmf17 L2/3.1','nmf78 L2/3.5','nmf27 L2/3.2','nmf45 L2/3',
                           'nmf84 L2/3','nmf62 HATA','nmf69 HATA',
                           'nmf83 VIP/HTR3A','nmf73 HTR3A','nmf50 VIP','nmf67 CXCL14','nmf35 LAMP5.CGE','nmf47 LAMP5.MGE',
                           'nmf88 C1QL1','nmf46 PV.FS','nmf74 SST','nmf60 SST/CORT'
                       )

plot<-plotExpression(sce.pats,features=c('nmf5 GC','nmf26 GC','nmf10 GC','nmf14 GC',
                                    'nmf52 MC','nmf11 CA3.1','nmf63 CA3.2','nmf61 CA2','nmf15 CA1',
                                    'nmf32 ProS','nmf40 Sub.1','nmf54 Sub.2','nmf22 L6.2','nmf65 L6/6b',
                                    'nmf53 L6b','nmf51 L5/6.1','nmf68 L5/6','nmf17 L2/3.1','nmf78 L2/3.5','nmf27 L2/3.2','nmf45 L2/3',
                                    'nmf84 L2/3','nmf62 HATA','nmf69 HATA',
                                    'nmf83 VIP/HTR3A','nmf73 HTR3A','nmf50 VIP','nmf67 CXCL14','nmf35 LAMP5.CGE','nmf47 LAMP5.MGE',
                                    'nmf88 C1QL1','nmf46 PV.FS','nmf74 SST','nmf60 SST/CORT'),
                     x="cellType", colour_by="cellType", point_alpha=1, point_size=.1,add_legend=F,ncol=4,scales='free_y')+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.3)+
    theme(axis.text.x = element_text(angle = 90))

pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','violins_ecs.pdf'),h=10,w=8)
rasterize(plot,layers = c('Point'),dpi=300)
dev.off()
