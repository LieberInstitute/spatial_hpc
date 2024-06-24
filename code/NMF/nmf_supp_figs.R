######snRNAseq mouse supp figures
library(SingleCellExperiment)
library(scater)
library(scran)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RcppML)
load('sce_subset_pats.rda')

####nmf pats vs annotation
data<-as.data.frame(colData(sce.subset2))
colnames(data)[match(c('nmf52','nmf11','nmf63','nmf61','nmf15',
                       'nmf32','nmf40','nmf54','nmf22','nmf65',
                       'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                       'nmf84','nmf62','nmf69'),colnames(data))]<-c(
                           'nmf52 MC','nmf11 CA3.1','nmf63 CA3.2','nmf61 CA2','nmf15 CA1',
                           'nmf32 ProS','nmf40 Sub.1','nmf54 Sub.2','nmf22 L6.2','nmf65 L6/6b',
                           'nmf53 L6b','nmf51 L5/6.1','nmf68 L5/6','nmf17 L2/3.1','nmf78 L2/3.5','nmf27 L2/3.2','nmf45 L2/3',
                           'nmf84 L2/3','nmf62 HATA','nmf69 HATA'
                       )

colnames(data)[match(c('nmf60','nmf74','nmf46','nmf88','nmf47',
                       'nmf35','nmf67','nmf50','nmf73','nmf83'),colnames(data))]<-c(
                           'nmf60 SST/CORT','nmf74 SST','nmf46 PV.FS','nmf88 C1QL1','nmf47 LAMP5.MGE',
                           'nmf35 LAMP5.CGE','nmf67 CXCL14','nmf50 VIP','nmf73 HTR3A','nmf83 VIP/HTR3A'
                       )
colnames(data)[match(c('nmf5','nmf10','nmf14','nmf26'),colnames(data))]<-c(
                           'nmf5 GC','nmf10 GC','nmf14 GC','nmf26 GC'
                       )



indices<-c(
    'nmf5 GC','nmf10 GC','nmf14 GC','nmf26 GC',
    'nmf52 MC','nmf11 CA3.1','nmf63 CA3.2','nmf61 CA2','nmf15 CA1',
    'nmf32 ProS','nmf40 Sub.1','nmf54 Sub.2','nmf22 L6.2','nmf65 L6/6b',
    'nmf53 L6b','nmf51 L5/6.1','nmf68 L5/6','nmf17 L2/3.1','nmf78 L2/3.5','nmf27 L2/3.2','nmf45 L2/3',
    'nmf84 L2/3','nmf62 HATA','nmf69 HATA',
    "nmf83 VIP/HTR3A","nmf73 HTR3A","nmf50 VIP","nmf67 CXCL14","nmf35 LAMP5.CGE",
    "nmf47 LAMP5.MGE", "nmf88 C1QL1", "nmf46 PV.FS", "nmf74 SST", "nmf60 SST/CORT"
)

data$annotation<-factor(data$annotation,levels=c("GC.1","GC.2",
                         "CA4","CA3.1","CA3.2","CA2",'CA1','PS.1','PS.2','Sub',
                         'L6.1','L6.2','L6.3','L6b',
                         'L5/Po.1','L5/Po.2',
                         'L2/3.1','L2/3.2',"L2/3.3","L2/3.4","L2/3.5","L2/3.6","L2/3.7",
                         "GABA.1",  "GABA.2",  "GABA.3",  "GABA.4",  "GABA.5"))
levels(data$annotation)[1]<-'GC'
levels(data$annotation)[2]<-'GC'



sce.pats<-SingleCellExperiment(assays=list(counts=t(data[,indices])))
colData(sce.pats)<-colData(sce.subset)
pdf(file=here::here('plots','figures','supp_figures','supp_fig_ecs_dotplot.pdf'),h=10.5,w=14)
create_custom_dot_plot(data, "annotation", indices, "", "NMF pattern",
                       "annotation", "proportion\nnuclei with\nnonzero\nweight",
                       "aggregate\npredicted\nnuclei-level\nweights")+
    theme(axis.text=element_text(size=12,color='black'),text=element_text(size=12,color='black'))
dev.off()





###supp fig sample-level factors
pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','nmf18_violin.pdf'),h=4,w=16)
p<-plotColData(sce,x='sample_id',y='nmf18',point_size=0.5,add_legend=T,colour_by='brnum')+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.3)+
    theme(axis.text.x = element_text(angle = 90,size=16),axis.text.y=element_text(size=16),
          text=element_text(size = 18))+facet_wrap(~sce$sort,scale='free_x',nrow=1)
ggrastr::rasterise(p,dpi=300,layers='Point')
dev.off()

pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','nmf18_superfine.pdf'),h=2,w=8.5)
p<-plotColData(sce,x='superfine.cell.class',y='nmf18',point_size=0.5,add_legend=F,colour_by='fine.cell.class')+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.3)+
    theme(axis.text.x = element_text(angle = 90,size=8),axis.text.y=element_text(size=8),
          text=element_text(size = 10))+scale_color_manual(values=sn.fine.palette)
ggrastr::rasterise(p,dpi=300,layers='Point')
dev.off()


###supp fig sample-level factors
pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','nmf16_violin.pdf'),h=4,w=16)
p<-plotColData(sce,x='sample_id',y='nmf16',point_size=0.5,add_legend=T,colour_by='brnum')+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.3)+
    theme(axis.text.x = element_text(angle = 90,size=16),axis.text.y=element_text(size=16),
          text=element_text(size = 18))+facet_wrap(~sce$sort,scale='free_x',nrow=1)
ggrastr::rasterise(p,dpi=300,layers='Point')
dev.off()

pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','nmf16_superfine.pdf'),h=2,w=8.5)
p<-plotColData(sce,x='superfine.cell.class',y='nmf16',point_size=0.5,add_legend=F,colour_by='fine.cell.class')+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.3)+
    theme(axis.text.x = element_text(angle = 90,size=8),axis.text.y=element_text(size=8),
          text=element_text(size = 10))+scale_color_manual(values=sn.fine.palette)
ggrastr::rasterise(p,dpi=300,layers='Point')
dev.off()

###supp fig sex factors
pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','nmf37_violin.pdf'),h=1.5,w=4.25)
p<-plotColData(sce,x='sex',y='nmf37',point_size=0.5,add_legend=T,colour_by='sex')+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.3)+
    theme(axis.text.x = element_text(angle = 90,size=8),axis.text.y=element_text(size=8),
          text=element_text(size = 10))
ggrastr::rasterise(p,dpi=300,layers='Point')
dev.off()

pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','nmf28_violin.pdf'),h=1.5,w=4.25)
p<-plotColData(sce,x='sex',y='nmf28',point_size=0.5,add_legend=T,colour_by='sex')+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.3)+
    theme(axis.text.x = element_text(angle = 90,size=8),axis.text.y=element_text(size=8),
          text=element_text(size = 10))
ggrastr::rasterise(p,dpi=300,layers='Point')
dev.off()

###tech variable heatmap
pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','heatmap_tech_vars.pdf'),h=2,w=8.4)
pheatmap(t(cor(t(x@h),as.matrix(colData(sce)[,c('sum','detected','subsets_Mito_percent')]))),
         breaks=seq(0,.786,.786/100),color = viridis::magma(100),
         fontsize = 10,cluster_rows=F)
dev.off()

###heatmap showing nmf94 across groups

pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','nmf94_byDonor.pdf'),h=2,w=4.25)
p<-plotColData(spe,x='brnum',y='nmf94',point_size=0.1,add_legend=T,colour_by='brnum')+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.3)+
    theme(axis.text.x = element_text(angle = 90,size=8),axis.text.y=element_text(size=8),
          text=element_text(size = 10))#+facet_wrap(~sce$choroid,scale='free_x',nrow=1)
ggrastr::rasterise(p,dpi=300,layers='Point')
dev.off()


pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','nmf94_byChoroid.pdf'),h=2.5,w=8.4)
p<-plotColData(sce,x='brnum',y='nmf94',point_size=0.1,add_legend=T,colour_by='choroid')+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.3)+
    theme(axis.text.x = element_text(angle = 90,size=8),axis.text.y=element_text(size=8),
          text=element_text(size = 10))+facet_wrap(~sce$choroid,scale='free_x',nrow=1)
ggrastr::rasterise(p,dpi=300,layers='Point')
dev.off()



####make barplot of # choroid cells by brnum
df <- as.data.frame(colData(sce))

# Filtering for choroid == TRUE
choroid_true_df <- df[df$choroid == TRUE, ]

# Counting the number of TRUE per brnum level
choroid_count <- table(choroid_true_df$brnum)

# Converting to dataframe for ggplot
choroid_count_df <- as.data.frame(choroid_count)
names(choroid_count_df) <- c("brnum", "count")


pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','choroid_barplot_by_donor.pdf'),h=2,w=4)
# Plotting
ggplot(choroid_count_df, aes(x = brnum, y = count,fill=brnum)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x = "donor", y = "number of choroid nuclei")+
    theme(axis.text.x = element_text(angle = 90,color='black'),
          axis.text.y = element_text(color='black'))
dev.off()



brains <- unique(spe$brnum)


speb <- spe[, (colData(spe)$brnum == brains[[6]])]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','nmf_94_spotplots.pdf'),h=7,w=7)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'nmf94',
    highlight = "broad.domain",
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    #palette = spatial.palette,
    image=F,
    values=spatial.palette3)+
    scale_fill_distiller(
        type = "seq",
        palette = rev('Greys'),
        direction=1)
dev.off()


pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','nmf_48_spotplots.pdf'),h=7,w=7)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'nmf30',
    highlight = "broad.domain",
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    #palette = spatial.palette,
    image=F,
    values=spatial.palette3)+
    scale_fill_distiller(
        type = "seq",
        palette = rev('Greys'),
        direction=1)
dev.off()


####multicolor plots for astro supp figure

speb <- spe[, (colData(spe)$brnum == brains[[3]])]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','astro_pattern_6423.pdf'),h=7,w=7)
plotVisiumRGB(speb,green='nmf76',blue='nmf81',pink='nmf79',image=F,highlight='neuron_cell_body',values=c('black','black'))
dev.off()

speb <- spe[, (colData(spe)$brnum == brains[[7]])]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','astro_pattern_8325.pdf'),h=8,w=7)
plotVisiumRGB(speb,green='nmf76',blue='nmf81',pink='nmf79',image=F,highlight='neuron_cell_body',values=c('black','black'))
dev.off()


###marker heatmap
pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','astro_heatmap.pdf'),h=2,w=4)
pheatmap(t(x@w[c('GFAP','SLC38A1','TNC','COL21A1','ADAMTSL3',
                 'SLC1A2','SLC1A3','SLC25A18','GPC5','NKAIN3',
                 'AQP4','APOE','ID4','AGT','GJA1'),
               c('nmf76','nmf79','nmf81')]),cluster_rows=F,cluster_cols=F,
         color=colorRampPalette(RColorBrewer::brewer.pal(9, "Greys"))(100),
         breaks=seq(0,0.0025,0.000025),scale='none',fontsize = 10)
dev.off()

###go dotplot astro
pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','astro_go_dotplot.pdf'),h=12,w=6)
merged<-list(go$nmf76,go$nmf79,go$nmf81)
names(merged)<-c('nmf76','nmf79','nmf81')
merged<-merge_result(merged)
dotplot(merged,showCategory=3)+scale_fill_distiller(
    type = "seq",
    palette = rev('Greys'),
    direction=-1)
dev.off()

###GO plots
go1<-go[c(2, 3, 4, 7, 8, 9, 10, 13, 14, 15, 17, 18, 19, 20, 21, 23, 24, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 64, 67, 68, 69, 70, 71, 72, 73, 76, 77, 78, 79, 80, 81, 82, 83, 84, 86, 88, 89, 90, 91, 92, 93)]
go2<-go[!names(go) %in% names(go1)]

go1<-merge_result(go1)
go2<-merge_result(go2)

pdf(file=here::here('plots','figures','supp_figures','supp_figures_nmf','noncelltype_goplots.pdf'),h=11,w=8.5)
dotplot(go2,showCategory=1)+theme(axis.text.y=element_text(size=6),axis.text.x=element_text(size=6,angle=90))
dev.off()



library(reshape2)
data<-as.data.frame(sce_pyr$superfine.cell.class)
colnames(data)<-'superfine.cell.class'
onehot_superfine.cell.class <-  dcast(data = data, rownames(data) ~ superfine.cell.class, length)
rownames(onehot_superfine.cell.class)<-onehot_superfine.cell.class[,1]
onehot_superfine.cell.class[,1]<-as.numeric(onehot_superfine.cell.class[,1])
onehot_superfine.cell.class<-onehot_superfine.cell.class[order(onehot_superfine.cell.class[,1],decreasing=F),]
onehot_superfine.cell.class[,1]<-NULL

weights<-as.data.frame(colData(sce_pyr)[,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                                           'nmf32','nmf40','nmf54','nmf22','nmf65',
                                           'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                                           'nmf84','nmf62','nmf69')])

#####Get all the projections in line!

###1: traced projections from Cembrowski et al:
colData(mouse)[,c(8:17)]<-NULL
mouse$projection<-factor(mouse$note2)
levels(mouse$projection)<-c('Amygdala','Nucleus accumbens','Prefrontal cortex')
##
# Translate from one species to the other using the orthology
names <- orthology[orthology$Column3 %in% rownames(mouse),]

names <- names[match(rownames(mouse), names$Column3),]

setdiff(names$Column3, rownames(mouse))

rownames(mouse) <- names$Column1


set.seed(1029)
i<-intersect(rownames(mouse),rownames(x@w))
loadings<-x@w
loadings<-loadings[rownames(loadings) %in% i,]
mouse2<-mouse[rownames(mouse) %in% i,]
loadings<-loadings[match(rownames(mouse2),rownames(loadings)),]
proj<-RcppML::project(loadings,logcounts(mouse2),L1=0)

col_sums<-colSums(t(proj))

# Rescale each column
proj<-t(proj)
proj2 <- apply(proj,1, function(row) row / col_sums)

# Check if columns now sum to 1
print(rowSums(proj2))

colData(mouse)<-cbind(colData(mouse),t(proj2))

##get aggregate projected weights
data<-as.data.frame(colData(mouse))

heat<-aggregate(data[,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                        'nmf32','nmf40','nmf54','nmf22','nmf65',
                        'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                        'nmf84','nmf62','nmf69')],
                by=list(data$projection),FUN=sum)
rownames(heat)<-heat$Group.1
heat<-heat[,-1]
pheatmap(heat,color=colorRampPalette(brewer.pal(n = 7, name =
                                                    "Greys"))(100),cluster_cols=F,cluster_rows=F)


###2: ding et al
sce.ding$subclass<-ifelse(sce.ding$cluster_id %in% c(1:3),'Sub-pyramidal',
                          ifelse(sce.ding$cluster_id %in% c(6:8),'ProS-pyramidal_deep',
                                 ifelse(sce.ding$cluster_id %in% c(10:14),'ProS-pyramidal_sup',
                                        ifelse(sce.ding$cluster_id %in% c(17:19),'HATA-pyramidal',
                                               ifelse(sce.ding$cluster_id %in% c(20:23),'ProS-pyramidal_mostSup',
                                                      ifelse(sce.ding$cluster_id %in% c(24:26),'Sub/ProS-polymorphic',
                                                             ifelse(sce.ding$cluster_id %in% c(28:29),'Sub/ProS/HATA-L6',
                                                                    ifelse(sce.ding$cluster_id %in% c(4:5),'ProS-Sub_border',
                                                                           ifelse(sce.ding$cluster_id %in% c(9,15,16),'CA1-pyramidal',
                                                                                  'PreS')))))))))
sce.ding$subclass<-factor(sce.ding$subclass)
sce.ding$subclass<-factor(sce.ding$subclass,
                          levels=levels(sce.ding$subclass)[c(1,7,6,4,
                                                             8,9,3,10,2,5)])
col_sums<-colSums(
    as.data.frame(colData(sce.ding)[,c(16:115)]))

# Rescale each column
colData(sce.ding)[,c(16:115)] <- t(apply(as.data.frame(colData(sce.ding)[,c(16:115)]),
                                         1, function(row) row / col_sums))

# Check if columns now sum to 1
print(colSums(as.data.frame(colData(sce.ding)[,c(16:115)])))

##get aggregate projected weights
data<-as.data.frame(colData(sce.ding))

heat<-aggregate(data[,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                        'nmf32','nmf40','nmf54','nmf22','nmf65',
                        'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                        'nmf84','nmf62','nmf69')],
                by=list(data$subclass),FUN=sum)
rownames(heat)<-heat$Group.1
heat<-heat[,-1]
pheatmap(heat,color=colorRampPalette(brewer.pal(n = 7, name =
                                                    "Greys"))(100),cluster_cols=F,cluster_rows=F)

Generally, seven major subclasses can be identified in the Sub-PS-HA region
(Figure 1C; Table S1): Spy (RC1–RC3), deep PSpy (RC6–RC8), superficial PSpy
(RC10–RC14), HApy (RC17–RC19), most superficial PSpy (RC20–RC23), Spo and PSpo
(RC24–RC26), and L6 of the Sub and HA (RC28 and RC29). Cells from two small
clusters (RC4 and RC5) lie in the PSpy at the border between the ventral Sub
and PS. Therefore, taking away adjoining CA1 (RC9, RC15, and RC16) and PrSd
(i.e., PoS; RC27) and adding L6b of the Sub and PS from the OC (OC48 and OC49),
a total of 27 clusters or cell types were revealed in the Sub-PS-HA region.

###3: Allen Institute
##rhistory stuff

##set up the anndata

##once you've got the keys set up:

##once you've got the onehots set up:
p<-cor(onehot_class,mat)
pheatmap(p,cluster_cols=F,cluster_rows=F)

p2<-heat[c('017 CA3 Glut','016 CA1-ProS Glut','022 L5 ET CTX Glut','023 SUB-ProS Glut',
           '032 L5 NP CTX Glut','033 NP SUB Glut','031 CT SUB Glut','030 L6 CT CTX Glut',
           '028 L6b/CT ENT Glut','005 L5 IT CTX Glut','006 L4/5 IT CTX Glut','004 L6 IT CTX Glut',
           '019 L2/3 IT PPP Glut','018 L2 IT PPP-APr Glut','008 L2/3 IT ENT Glut','007 L2/3 IT CTX Glut',
           '011 L2 IT ENT-po Glut','014 LA-BLA-BMA-PA Glut','012 MEA Slc17a7 Glut'),]


##get aggregate projected weights
data<-as.data.frame(colData(mch2))

heat<-aggregate(data[,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                        'nmf32','nmf40','nmf54','nmf22','nmf65',
                        'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                        'nmf84','nmf62','nmf69')],
                by=list(data$Target),FUN=sum)
rownames(heat)<-heat$Group.1
heat<-heat[,-1]
pheatmap(heat,color=colorRampPalette(brewer.pal(n = 7, name =
                                                    "Greys"))(100),cluster_cols=F,cluster_rows=F)


##########umaps#############
sce<-makeVisiumRGB(sce,vars=c('nmf11','nmf15','nmf52','nmf61'))

pdf(file=here::here('plots','figures','figure_6','UMAP_rgb_hpc.pdf'),h=5,w=5)
p<-plotUMAP(sce,colour_by='RGB',text_by='fine.cell.class',point_size=0.01)+scale_color_identity()
ggrastr::rasterize(p,dpi=500,layers='Point')
dev.off()

sce<-makeVisiumRGB(sce,vars=c('nmf54','nmf17','nmf32','nmf40'))

pdf(file=here::here('plots','figures','figure_6','UMAP_rgb_sub.pdf'),h=5,w=5)
p<-plotUMAP(sce,colour_by='RGB',text_by='fine.cell.class',point_size=0.01)+scale_color_identity()
ggrastr::rasterize(p,dpi=500,layers='Point')
dev.off()



###marker gene heatmap
heat<-x@w[rownames(x@w) %in% pick,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                                    'nmf32','nmf40','nmf54','nmf22','nmf65',
                                    'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                                    'nmf84','nmf62','nmf69')]

pheatmap(heat,color=colorRampPalette(brewer.pal(n = 7, name =
                                                    "Greys"))(100),cluster_cols=F,cluster_rows=T)



glut<-c('GRIK1','GRIK3','GRIK4','GRIA4','GRID2','GRM8','GRM3')
gad<-c('GABRB2','GABRA2','GABRA1')
htr<-c('HTR2C','HTR4','HTR7','HTR1E','HTR2A')

nt<-c(glut,gad,htr)

genes<-genes[names(genes) %in% c('nmf52','nmf11','nmf63','nmf61','nmf15',
                                 'nmf32','nmf40','nmf54','nmf22','nmf65',
                                 'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                                 'nmf84','nmf62','nmf69')]



#genes in diff cats:
ion<-c('')
ntr<-
    np<-

    ####get pyr-spec markers
    pyr<-c('nmf52','nmf11','nmf63','nmf61','nmf15',
           'nmf32','nmf40','nmf54','nmf22','nmf65',
           'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
           'nmf84','nmf62','nmf69')
###get markers for ORA after filtering mito genes and non-protein coding genes
loads<-x@w
loads<-loads[,colnames(loads) %in% pyr]
##get rid of genes with no weights in any factor
no_expr <- which(rowSums(loads) == 0)
# length(no_expr)
# length(no_expr) / nrow(loads) * 100
loads <- loads[-no_expr, ]
##now filter mito genes and non-protein coding genes
protein<-rownames(sce)[rowData(sce)$gene_type=='protein_coding']
loads<-loads[rownames(loads) %in% protein,]
mito<-rownames(sce)[which(seqnames(sce) == "chrM")]
loads<-loads[!rownames(loads) %in% mito,]
##now get markers
marks<-patternMarkers(loads,x@h,'all',1,100)



markers<-marks$PatternMarkers[c(52,11,63,61,15,32,40,54,22,65,53,
                                51,68,17,78,27,45,84,62,69)]



features<-c('CARTPT','ZNF208','HGF','TGFB2','FIBCD1','TG','FN1',
            'COL12A1','SULF1','SDK2','GSG1L','RORB','AGBL1','VIPR2','NIPAL2',
            'CRHR2','CBLN4','RSPO2','KLHL1','TACR3')



heat<-x@w[features,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                     'nmf32','nmf40','nmf54','nmf22','nmf65',
                     'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                     'nmf84','nmf62','nmf69')]

pdf(file=here::here('plots','figures','figure_6','marker_gene_heatmap.pdf'),h=4.5,w=4)
pheatmap(heat,cluster_cols=F,cluster_rows=F,
         color=colorRampPalette(brewer.pal(n = 7, name ="Greys"))(100),
         breaks=seq(0,0.001,0.00001
         ))
dev.off()

library(dplyr)
library(tidyr)
library(ggplot2)

# Sample data (replace this with your actual data frame)
# data <- your_data_frame
data<-as.data.frame(colData(mch2))
data<-as.tibble(data)

# Specified order for the x-axis
features_order <- c('nmf52','nmf11','nmf63','nmf61','nmf15',
                    'nmf32','nmf40','nmf54','nmf22','nmf65',
                    'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                    'nmf84','nmf62','nmf69')

# Selecting relevant columns and melting the data into long format
long_data <- data %>%
    dplyr::select(Subclass, all_of(features_order)) %>%
    pivot_longer(cols = -Subclass, names_to = "Feature", values_to = "Value")

# Calculating sum and proportion of non-zero cells
stats <- long_data %>%
    group_by(Subclass, Feature) %>%
    summarize(
        Sum = sum(Value),
        Mean=mean(Value),
        NonZeroProportion = sum(Value != 0) / n(),
        total=sum(Value != 0)  # Explicit proportion calculation
    ) %>%
    ungroup()

# Adjusting Feature as a factor with the specified order
stats$Feature <- factor(stats$Feature, levels = features_order)

# Creating the plot with switched axes
ggplot(stats, aes(x = Feature, y = Subclass, size = NonZeroProportion, color = Sum)) +
    geom_point() +
    scale_size_continuous(range = c(0, 12)) + # Adjust point sizes as needed
    scale_color_gradient(low = "white", high = "black") + # Greyscale color scale
    theme_minimal() +
    labs(
        title = "Modified Dot Plot",
        x = "Feature",
        y = "Subclass",
        size = "Proportion of Non-Zero Cells",
        color = "Sum of Values"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust x-axis text angle for readability


###Target
library(dplyr)
library(tidyr)
library(ggplot2)

# Sample data (replace this with your actual data frame)
# data <- your_data_frame

# Specified order for the x-axis
features_order <- c('nmf52','nmf11','nmf63','nmf61','nmf15',
                    'nmf32','nmf40','nmf54','nmf22','nmf65',
                    'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                    'nmf84','nmf62','nmf69')

# Selecting relevant columns and melting the data into long format
long_data <- data %>%
    dplyr::select(Target, all_of(features_order)) %>%
    pivot_longer(cols = -Target, names_to = "Feature", values_to = "Value")

# Calculating sum and proportion of non-zero cells
stats <- long_data %>%
    group_by(Target, Feature) %>%
    summarize(
        Sum = sum(Value),
        Mean=mean(Value),
        NonZeroProportion = sum(Value != 0) / n(),
        total=sum(Value != 0)
    ) %>%
    ungroup()

# Adjusting Feature as a factor with the specified order
stats$Feature <- factor(stats$Feature, levels = features_order)

# Creating the plot with switched axes
ggplot(stats, aes(x = Feature, y = Target, size = Sum, color = NonZeroProportion)) +
    geom_point() +
    scale_size_continuous(range = c(0, 12)) + # Adjust point sizes as needed
    scale_color_gradient(low = "white", high = "black") + # Greyscale color scale
    theme_minimal() +
    labs(
        title = "Modified Dot Plot",
        x = "Feature",
        y = "Target",
        size = "Proportion of Non-Zero Cells",
        color = "Sum of Values"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust x-axis text angle for readability



####supp figure aggregate weights heatmap
##get aggregate projected weights
data<-as.data.frame(colData(sce))

heat<-aggregate(data[,paste0('nmf',c(1:100))],
                by=list(sce$round),FUN=sum)
rownames(heat)<-heat$Group.1
heat<-heat[,-1]

pdf(file=here::here('plots','figures','supp_figures',
                    'supp_figures_nmf','superfine_heatmap.pdf'),h=16,w=20)
pheatmap(heat,
         color=colorRampPalette(brewer.pal(n = 7, name =
                                               "Greys"))(100),cluster_cols=T,cluster_rows=F,
         clustering_distance_cols = 'manhattan',legend=T)

dev.off()
###supp figure round correlation

pdf(file=here::here('plots','figures','supp_figures',
                    'supp_figures_nmf','round_heatmap.pdf'),h=6,w=17)
pheatmap(cor(onehot_round,t(x@h)),cluster_rows=F,cluster_cols=F)

dev.off()

###supp figure round correlation

pdf(file=here::here('plots','figures','supp_figures',
                    'supp_figures_nmf','sex_heatmap.pdf'),h=4,w=17)
pheatmap(cor(onehot_sex,t(x@h)),cluster_rows=F,cluster_cols=F)

dev.off()

###supp figure round correlation

pdf(file=here::here('plots','figures','supp_figures',
                    'supp_figures_nmf','sex_heatmap.pdf'),h=4,w=17)
pheatmap(cor(onehot_brnum,t(x@h)),cluster_rows=F,cluster_cols=F)

dev.off()



##get aggregate projected weights
data<-as.data.frame(colData(sce))

heat<-aggregate(data[,paste0('nmf',c(1:100))],
                by=list(sce$superfine.cell.class),FUN=mean)
rownames(heat)<-heat$Group.1
heat<-heat[,-1]
pheatmap(heat,color=colorRampPalette(brewer.pal(n = 7, name =
                                                    "Greys"))(100),cluster_cols=T,cluster_rows=T)


###for supp figure--markers vs nmf proj
sce<-sce[,!sce$fine.cell.class %in% c('Amy','Thal','GABA.PENK','GC','Cajal')]
sce<-sce[,sce$mid.cell.class=='ExcN']
sce<-sce[,!sce$fine.cell.class=='Cajal']

sce$superfine.cell.class<-droplevels(sce$superfine.cell.class)

sce<-sce[rownames(sce) %in% rownames(mch2),]
mch2<-mch2

markers<-get_mean_ratio2(sce,cellType_col = 'superfine.cell.class',
                         assay_name='logcounts')

library(dplyr)

df<-markers
# Filter the top 10 genes by rank_ratio for each cellType.target
topGenes <- df %>%
    arrange(cellType.target, rank_ratio) %>%
    group_by(cellType.target) %>%
    slice_head(n = 5) %>%
    ungroup()

# Split the data frame into a list of data frames, one for each cellType.target
df_list <- split(topGenes,
                 topGenes$cellType.target)

print(df_list)

sumGenes<-function(x){
    if(nrow(x)==1){
        summed<-assay(mch2,'X')[x$gene,]
        return(summed)
    }
    else{
        print('summing!')
        summed<-colSums(assay(mch2,'X')[x$gene,])
        return(summed)
    }}

sums<-lapply(df_list,sumGenes)
sums<-as.data.frame(sums)
identical(rownames(sums),colnames(mch2))
###TRUE
# colData(mch2)<-colData(mch2)[,!colnames(colData(mch2)) %in% indices]

#####cbind
colData(mch2)<-cbind(colData(mch2),sums)

###make dotplot
data<-as.data.frame(colData(mch2))

heat<-aggregate(data[,indices],
                by=list(data$Subclass),FUN=sum)
rownames(heat)<-heat$Group.1
heat<-heat[,-1]

pdf(file=here::here('plots','figures','supp_figures',
                    'supp_figures_nmf','epiretro_5markers.pdf'),h=3,w=4)
pheatmap(heat,cluster_cols=F,cluster_rows=F,scale='row',fontsize = 8)
dev.off()





heat<-aggregate(data[,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                        'nmf32','nmf40','nmf54','nmf22','nmf65',
                        'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                        'nmf84','nmf62','nmf69')],
                by=list(data$Subclass),FUN=sum)
rownames(heat)<-heat$Group.1
heat<-heat[,-1]
pheatmap(heat,color=colorRampPalette(brewer.pal(n = 7, name =
                                                    "Greys"))(100),cluster_cols=F,cluster_rows=F)

