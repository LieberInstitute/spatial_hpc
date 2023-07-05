data<-as.data.frame(spez$ManualAnnotation)
colnames(data)<-'ManualAnnotation'
onehot_manualannotation <-  dcast(data = data, rownames(data) ~ ManualAnnotation, length)
rownames(onehot_manualannotation)<-onehot_manualannotation[,1]
onehot_manualannotation[,1]<-as.numeric(onehot_manualannotation[,1])
onehot_manualannotation<-onehot_manualannotation[order(onehot_manualannotation[,1],decreasing=F),]
onehot_manualannotation[,1]<-NULL

plotGroupedHeatmap(spe,group = spe$PRECAST_k17, features='SLC17A7')

rownames(spe)<-
    uniquifyFeatureNames(rowData(spe)$gene_id, rowData(spe)$gene_name)


brains <- unique(spe$brnum)


speb <- spe[, (colData(spe)$brnum == brains[2])]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

#pdf(here::here('plots','figure_1','visium_FIBCD1.pdf'),h=6,w=6)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'EPHA1',
    highlight = NULL,
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    #palette = as.vector(alphabet.colors()),
    image=T
)+scale_fill_viridis(option='inferno')


pdf(here::here('plots','figure_1','visium_histology.pdf'),h=4,w=4)
plotVisium(
    speb,
    spots = F,
    fill = NULL,
    highlight = NULL,
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    # palette = c('blue','red'),
    image=T)

plotExpression(spe,features="SCG2",
               x="cluster", colour_by="cluster", point_alpha=0.5, point_size=.7,add_legend=F)+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.3)+
    theme(axis.text.x = element_text(angle = 45),
          text=element_text(size = 13))+
    labs(x='Cell Type:Condition',y='log2 normalized counts')

spe$cluster<-as.character(spe$PRECAST_k17_nnSVG_rep)
spe$cluster<-ifelse(spe$cluster %in% 17,'CP',
                    ifelse(spe$cluster %in% 12,'Vasc',
                           ifelse(spe$cluster %in% 4,'Olig.1',
                                  ifelse(spe$cluster %in% 11,'Olig.2',
                                         ifelse(spe$cluster %in% 1,'Olig.3',
                                                ifelse(spe$cluster %in% 9,'SLM/SR',
                                                       ifelse(spe$cluster %in% 6,'SR/SL',
                                                              ifelse(spe$cluster %in% 13,'DG-ML',
                                                                     ifelse(spe$cluster %in% 5,'SGZ/SO/Sub-ML',
                                                                            ifelse(spe$cluster %in% 3,'GABA',
                                                                                   ifelse(spe$cluster %in% 8,'GCL',
                                                                                          ifelse(spe$cluster %in% 7,'CA4/CA3',
                                                                                                 ifelse(spe$cluster %in% 15,'CA3',
                                                                                                        ifelse(spe$cluster %in% 10,'CA2/CA1',
                                                                                                               ifelse(spe$cluster %in% 14,'CA1',
                                                                                                                      ifelse(spe$cluster %in% 2,'Sub',
                                                                                                                             'EC'))))))))))))))))
plotGroupedHeatmap(spe,group='PRECAST_k17_nnSVG_rep',features=c('SNAP25','NRGN',
                                                                'MT3','MT-CO1',
                                                                'TTR','CLDN5',
                                                                'MBP','PLP1'))
spe$cluster<-factor(spe$cluster,levels=c('CP','Vasc','Olig.1','Olig.2','Olig.3',
                                         'SLM/SR','SR/SL','DG-ML','SGZ/SO/Sub-ML',
                                         'GABA','GCL','CA4/CA3','CA3','CA2/CA1','CA1','Sub','EC'))
features=c('PRLR','TAGLN','MOBP','OLIG1','APOC1','ETNPPL','SYT1','SST','GAD1',
           'PPFIA2','CARTPT','KCNQ5','FIBCD1','COL5A2','ETV1','MEF2C')
plotDots(spe,group='cluster',
         features=rev(features),
         colour=c('white','black'))+theme(axis.text.x=element_text(angle=45))

plotDots(spe,group='PRECAST_k17_nnSVG_rep',
         features=c('FOLR1','TAGLN','FIBCD1','GAD2',
                    'MOBP','ETV1','SST','APOC1','CARTPT','TRHDE','OSR1'),
         colour=c('white','red'))



set.seed(1325)
sce3<-runUMAP(sce,dimred='PCA')

tab$type<-rep(NA,51)
tab$type[tab$annoType %in% c('Amy.1','Amy.2','Amy.3')]<-'Amygdala'
tab$type[tab$annoType %in% c("Astro.1","Astro.2")]<-'Astro'
tab$type[tab$annoType %in% c("CA1","CA2","CA3",'CA3.deep','MC')]<-'CA1-4'
tab$type[tab$annoType %in% c("CP.1",'CP.2')]<-'Choroid'
tab$type[tab$annoType %in% c("Ependy")]<-'Ependymal'
tab$type[tab$annoType %in% c("HATA")]<-'HATA'
tab$type[tab$annoType %in% paste0("GC.",c(1:4))]<-'GC'
tab$type[tab$annoType %in% c("Endo/PC",'SMC','VLMC')]<-'Vascular'
tab$type[tab$annoType %in% c("OPC",'COP')]<-'OPC'
tab$type[tab$annoType %in% c("Oligo.1",'Oligo.2','Oligo.NB')]<-'Oligo'
tab$type[tab$annoType %in% c("Micro",'Tcell/Macro')]<-'Immune'
tab$type[tab$annoType %in% c("RHP.L2.1", "RHP.L2.2", "RHP.L2.3", "RHP.L2.4",
                             "RHP.L2.5", "RHP.L3", "RHP.L5", "RHP.L6.1", "RHP.L6.2",
                             "RHP.L6b")]<-'RHP'
tab$type[tab$annoType %in% c("Sub.Deep",'Sub.Sup')]<-'Sub'
tab$type[tab$annoType %in% c("Thal")]<-'Thal'
tab$type[tab$annoType %in% c("Cajal")]<-'Cajal-Retzius'
tab$type[tab$annoType %in% c("CORT", "HTR3A","LAMP5 KIT", "LAMP5 LHX6","MEIS2",
                             "NKX2-1", "PVALB", "SST",
                             "VIP")]<-'GABA'

set.seed(1029)
sce2<-runUMAP(sce,dimred='HARMONY')



71: >50
55: >50
61: >50
9: >50
21: >50
35: >100
49: >70


spe$type71<-ifelse(spe$nmf71>=50,1,0)
spe$type55<-ifelse(spe$nmf55>=50,1,0)
spe$type61<-ifelse(spe$nmf61>=50,1,0)
spe$type9<-ifelse(spe$nmf9>=50,1,0)
spe$type21<-ifelse(spe$nmf21>=30,1,0)
spe$type35<-ifelse(spe$nmf35>=75,1,0)
spe$type49<-ifelse(spe$nmf49>=75,1,0)

brains<-unique(spe$brnum)


speb <- spe[, (colData(spe)$brnum == brains[2])]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

plotVisium(
    speb,
    spots = TRUE,
    fill = 'ENHO',
    highlight = NULL,
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    # palette = c('blue','red'),
    image=T)+scale_fill_viridis(option='inferno')

cols<-cols[-c(50,73,24,26,40)]

tab<-data.frame('cluster'=1:18,'annotation2'=rep(NA,18),'annotation'=rep(NA,18))
tab$annotation[tab$cluster %in% c(4)]<-'GCL'
tab$annotation[tab$cluster %in% c(1)]<-'SUB.RHP.2'
tab$annotation[tab$cluster %in% c(2)]<-'DG.ML'
tab$annotation[tab$cluster %in% c(3)]<-'CA1.1'
tab$annotation[tab$cluster %in% c(5)]<-'SR.SLM'
tab$annotation[tab$cluster %in% c(6)]<-'WM.1'
tab$annotation[tab$cluster %in% c(7)]<-'CA2.4.1'
tab$annotation[tab$cluster %in% c(8)]<-'RHP.1'
tab$annotation[tab$cluster %in% c(9)]<-'CA1.2'
tab$annotation[tab$cluster %in% c(10)]<-'GABA'
tab$annotation[tab$cluster %in% c(11)]<-'Vascular'
tab$annotation[tab$cluster %in% c(12)]<-'Choroid'
tab$annotation[tab$cluster %in% c(13)]<-'SL.SR'
tab$annotation[tab$cluster %in% c(14)]<-'SUB.RHP.1'
tab$annotation[tab$cluster %in% c(15)]<-'WM.2'
tab$annotation[tab$cluster %in% c(16)]<-'WM.3'
tab$annotation[tab$cluster %in% c(17)]<-'CA2.4.2'
tab$annotation[tab$cluster %in% c(18)]<-'WM.4'


tab$annotation2[tab$cluster %in% c(4)]<-'Neuron'
tab$annotation2[tab$cluster %in% c(1)]<-'Neuron'
tab$annotation2[tab$cluster %in% c(2)]<-'Neuropil'
tab$annotation2[tab$cluster %in% c(3)]<-'Neuron'
tab$annotation2[tab$cluster %in% c(5)]<-'Neuropil'
tab$annotation2[tab$cluster %in% c(6)]<-'WM'
tab$annotation2[tab$cluster %in% c(7)]<-'Neuron'
tab$annotation2[tab$cluster %in% c(8)]<-'Neuron'
tab$annotation2[tab$cluster %in% c(9)]<-'Neuron'
tab$annotation2[tab$cluster %in% c(10)]<-'Neuron'
tab$annotation2[tab$cluster %in% c(11)]<-'Vascular'
tab$annotation2[tab$cluster %in% c(12)]<-'Vascular'
tab$annotation2[tab$cluster %in% c(13)]<-'Neuropil'
tab$annotation2[tab$cluster %in% c(14)]<-'Neuron'
tab$annotation2[tab$cluster %in% c(15)]<-'WM'
tab$annotation2[tab$cluster %in% c(16)]<-'WM'
tab$annotation2[tab$cluster %in% c(17)]<-'Neuron'
tab$annotation2[tab$cluster %in% c(18)]<-'WM'

spe$cluster<-tab$annotation[match(spe$nnSVG_PRECAST_k18,tab$cluster)]
spe$broad<-factor(tab$annotation2[match(spe$nnSVG_PRECAST_k18,tab$cluster)])


spe$cluster2<-factor(spe$cluster,levels=c("GCL","CA2.4.1", "CA2.4.2",
                                          "CA1.1", "CA1.2","SUB.RHP.1", "SUB.RHP.2",
                                          "RHP.1","GABA","SL.SR","DG.ML", "SR.SLM",
                                          paste0("WM.",c(2,1,3,4)),"Vascular","Choroid"))


brains <- unique(spe$brnum)


speb <- spe[, (colData(spe)$brnum == brains[2])]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

pdf(here::here('plots','figures','figure_2','visium_clusters.pdf'),h=8,w=8)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'S1PR1',
    highlight = NULL,
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    #palette = as.vector(palette36.colors(18)),
    image=T
)
dev.off()


features<-c('GRP','KCNQ5','FNDC1','POU3F1','GFRA1','AKAP12','KRT17','MEF2C','SLC32A1',
            'APOC1','MOBP','BCAS1','SFRP2','ACTA2','PRLR')

'ETNPPL','MOBP','GPR17','PDGFRA','C3','SKAP1','CFAP73','VWF',
'CEMIP','ACTA2','FOLR1','SLC17A7','TNNT2','ABI3BP','TG',
'ST18','RGS14','FIBCD1','FN1','LPO','TLE4','CCN2','RORB','PCP4','CUX2','ESR1','NDNF',
'SHOX2','GAD2','LHX6','SST','CORT','BTBD11','CRABP1','LAMP5','CHST9','KIT','VIP','HTR3A','MEIS2')


pdf(file=here::here('plots','figures','figure_2','visium_dotplot.pdf'),h=5,w=7)
plotDots(spe,group='cluster',features=features,color=c('white','black'))+
    scale_y_discrete(limits=rev(features)) +
    #scale_x_discrete(limits=rev(levels(spe$cluster)))+
    theme(axis.text.x = element_text(angle = 90,vjust=0.75),
          text = element_text(size = 14))+labs(x= "Cluster", y = "Gene")
dev.off()


plotGroupedHeatmapMod <- function(spe, group, features, assay.type = "logcounts", color = pheatmap::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(25), show_colnames = F) {
    # Check inputs
    if (!group %in% colnames(colData(spe))) {
        stop(paste("Group", group, "not found in colData of the SingleCellExperiment object."))
    }

    if (!all(features %in% rownames(spe))) {
        stop(paste("Some or all features not found in the rownames of the SingleCellExperiment object."))
    }

    if (!assay.type %in% names(assays(spe))) {
        stop(paste("Assay type", assay.type, "not found in the assays of the SingleCellExperiment object."))
    }

    # Subset spe to features of interest
    spe_subset <- spe[features, ]

    # Aggregate assay values (e.g. logcounts) across cells (by group)
    avg_logcounts <- aggregateAcrossCells(spe_subset, DataFrame(
        cluster = colData(spe_subset)[['cluster']],
        sample_id = spe$sample_id))

    # Convert to matrix
    avg_logcounts <- as.matrix(assays(avg_logcounts)[[1]])

    # Remove potential NA values
    avg_logcounts[is.na(avg_logcounts)] <- 0

    # Prepare annotation data
    annotation_data <- data.frame(broad = colnames(avg_logcounts))
    rownames(annotation_data) <- colnames(avg_logcounts) # The rownames of annotation_data should match the column names of avg_logcounts

    # Plot
    pheatmap::pheatmap(avg_logcounts, color = color, show_colnames = show_colnames, annotation_col = annotation_data)
}


pheatmap::pheatmap(avg_logcounts, color = color, show_colnames = T,scale='row')




ggplot(df_summary, aes(x = cluster, y = median_expression, color = gene, group = gene)) +
    geom_line() +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.line = element_line(color = "black"),
          legend.position = "none",
          axis.text.x = element_text(angle = 90),
          text = element_text(size = 10)) +
    labs(y = "Median Expression") +
    scale_color_viridis_d()


features<-c('ETNPPL','MOBP','GPR17','PDGFRA','C3','SKAP1','CFAP73','VWF',
            'CEMIP','ACTA2','FOLR1','SLC17A7','SLC17A6','TNNT2','POSTN','NPAS4','ABI3BP','TG',
            'ST18','RGS14','FIBCD1','FN1','LPO','TLE4','IGFBP3','CCN2','RORB','TRABD2A','PCP4','CUX2','QRFPR',
            'CRHR2','CBLN2','SNHG31','SH3RF2','MET','ESR1','DCSTAMP','NPFFR2','PAPPA2','NDNF',
            'SHOX2','GAD2','LHX6','SST','CORT','BTBD11','CRABP1','LAMP5','KIT','VIP','HTR3A','MEIS2')

features<-c('ETNPPL','MOBP','GPR17','PDGFRA','C3','SKAP1','CFAP73','VWF',
            'CEMIP','ACTA2','FOLR1','SLC17A7','TNNT2','ABI3BP','TG',
            'ST18','RGS14','FIBCD1','FN1','LPO','TLE4','CCN2','RORB','PCP4','CUX2','ESR1','NDNF',
            'SHOX2','GAD2','LHX6','SST','CORT','BTBD11','CRABP1','LAMP5','CHST9','KIT','VIP','HTR3A','MEIS2')

features=rownames(marks2[[2]][26:50,])
plotDots(spe,group='cluster',features=features,color=c('white','black')) +
    #scale_y_discrete(limits=rev(features)) +
    scale_x_discrete(limits=rev(levels(spe$cluster)))+
    theme(axis.text.x = element_text(angle = 90,vjust=0.75),
          text = element_text(size = 14))+ coord_flip()+
            labs(x= "Cell type", y = "Gene")
features=c('THY1','SYT1','CHN1','PLP1','GFAP','VIM'
#pdf(file=here::here('plots','figures','figure_2','violin_plots.pdf'),h=3,w=4.3)
#features=c('MDGA1','PRLR','MN1','PAPPA2','CDH23','NTS')
plotExpression(spe,features=c('SYT1','PLP1','AQP4','VIM'),
               x="cluster", colour_by="broad", point_alpha=0.5,add_legend=T,ncol=1)+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.3)+
    theme(axis.text.x = element_text(angle = 90),
          text=element_text(size = 11),
        axis.text=element_text(size=9))+theme(strip.text = element_blank())

    #scale_color_viridis_d()+
    labs(y='log2 normalized counts',x='Cluster')
#dev.off()

plotColData(spe,x='cluster',y='subsets_Mito_percent',colour_by='cluster',point_size=0.5,add_legend=F)+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.3)+
    theme(axis.text.x = element_text(angle = 90),
          text=element_text(size = 13))

ggcells(spe, mapping=aes(x=cluster, y=VIM,fill=broad)) +
    geom_boxplot()+theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'))+
    theme(legend.position='bottom',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
x2<-ggcells(spe, mapping=aes(x=cluster, y=PLP1,fill=broad)) +
    geom_boxplot()+theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
x3<-ggcells(spe, mapping=aes(x=cluster, y=VIM,fill=broad)) +
    geom_boxplot()+theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 13,colour='black'))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))

# First, prepare data for plotting
df <- assays(spe)$logcounts[genes, ]
df <- as.data.frame(t(df))
df$cluster <- colData(spe)$cluster
df <- df %>% gather(key = "gene", value = "expression", -cluster)

# Compute median expression in each cluster for each gene
df_summary <- df %>%
    group_by(cluster, gene) %>%
    summarise(median_expression = median(expression, na.rm = TRUE))

# Plot
ggplot(df_summary, aes(x = cluster, y = median_expression, color = gene, group = gene)) +
    geom_line() +
    theme(axis.text.x = element_text(angle = 90),
          text=element_text(size = 10)) +
    labs(y = "Median Expression") +
    scale_color_viridis_d()



c('PDGFRA','GPR17','SOX10','OLIG1','OLIG2','MAG','PLP1','CNP','MBP','MOBP','OPALIN','MOG')


#set.seed(602)
#message("running TSNE - ", Sys.time())
#sce <- runTSNE(sce, dimred = "HARMONY")
reducedDim(sce,'MNN')<-mnn$corrected
set.seed(0719)
message("running UMAP - ", Sys.time())
sce2 <- runUMAP(sce, dimred = "HARMONY")
message("running buildSNNGraph - ", Sys.time())
snn.gr <- buildSNNGraph(sce, k = 5, use.dimred = "HARMONY",type='jaccard')


message("running louvain - ", Sys.time())
set.seed(100)
clust_5 <- igraph::cluster_louvain(snn.gr,resolution=1.2)$membership
table(clust_5)
sce$cluster_harmony<-factor(clust_5)



tab<-data.frame('cluster'=1:51,'annotation'=rep(NA,51),'annoType'=rep(NA,51))
########Neurons
##Granule cells
tab$annotation[tab$cluster %in% c(10,28,39,40)]<-'EXC'
tab$annoType[tab$cluster %in% c(10,28,39,40)]<-'GC'
##Mossy cells
tab$annotation[tab$cluster %in% c(18)]<-'EXC'
tab$annoType[tab$cluster %in% c(18)]<-'MC'
##CA3
tab$annotation[tab$cluster %in% c(9)]<-'EXC'
tab$annoType[tab$cluster %in% c(9)]<-'CA3_1'
##CA3
tab$annotation[tab$cluster %in% c(33)]<-'EXC'
tab$annoType[tab$cluster %in% c(33)]<-'CA3_2'
##CA2
tab$annotation[tab$cluster %in% c(36)]<-'EXC'
tab$annoType[tab$cluster %in% c(36)]<-'CA2'
##CA1
tab$annotation[tab$cluster %in% c(11)]<-'EXC'
tab$annoType[tab$cluster %in% c(11)]<-'CA1'

##Sub_sup
tab$annotation[tab$cluster %in% c(23)]<-'EXC'
tab$annoType[tab$cluster %in% c(23)]<-'Sub.Sup'
##Sub_deep
tab$annotation[tab$cluster %in% c(19)]<-'EXC'
tab$annoType[tab$cluster %in% c(19)]<-'Sub.Deep'

##RHP
#tab$annotation[tab$cluster %in% c(25,29,38,43,45)]<-'EXC'
#tab$annoType[tab$cluster %in% c(25,29,38,43,45)]<-paste0('RHP.L2.',c(1:5))

tab$annotation[tab$cluster %in% c(5)]<-'EXC'
tab$annoType[tab$cluster %in% c(5)]<-'RHP.L2.3_1'

tab$annotation[tab$cluster %in% c(25,29,38,43,45)]<-'EXC'
tab$annoType[tab$cluster %in% c(25,29,38,43,45)]<-paste0('RHP.L2.3_',c(2:6))

tab$annotation[tab$cluster %in% c(21)]<-'EXC'
tab$annoType[tab$cluster %in% c(21)]<-'RHP.L5'

tab$annotation[tab$cluster %in% c(12,14)]<-'EXC'
tab$annoType[tab$cluster %in% c(12,14)]<-paste0('RHP.L6_',c(1,2))

tab$annotation[tab$cluster %in% c(20)]<-'EXC'
tab$annoType[tab$cluster %in% c(20)]<-'RHP.L6b'

tab$annotation[tab$cluster %in% c(4)]<-'EXC'
tab$annoType[tab$cluster %in% c(4)]<-'Thal'



##Amygdala
tab$annotation[tab$cluster %in% c(41,42,48)]<-'EXC'
tab$annoType[tab$cluster %in% c(41,42,48)]<-'Amy'
##HATA
tab$annotation[tab$cluster %in% c(32)]<-'EXC'
tab$annoType[tab$cluster %in% c(32)]<-'HATA'

##Cajal Retzius
tab$annotation[tab$cluster %in% c(46)]<-'EXC'
tab$annoType[tab$cluster %in% c(46)]<-'Cajal'

##GABAergics
tab$annotation[tab$cluster %in% c(47)]<-'INH'
tab$annoType[tab$cluster %in% c(47)]<-'MEIS2'

tab$annotation[tab$cluster %in% c(13)]<-'INH'
tab$annoType[tab$cluster %in% c(13)]<-'PV'

tab$annotation[tab$cluster %in% c(8)]<-'INH'
tab$annoType[tab$cluster %in% c(8)]<-'SST'

tab$annotation[tab$cluster %in% c(49)]<-'INH'
tab$annoType[tab$cluster %in% c(49)]<-'CRABP1'

tab$annotation[tab$cluster %in% c(3)]<-'INH'
tab$annoType[tab$cluster %in% c(3)]<-'LAMP5.LHX6'

tab$annotation[tab$cluster %in% c(15)]<-'INH'
tab$annoType[tab$cluster %in% c(15)]<-'LAMP5.KIT'

tab$annotation[tab$cluster %in% c(34)]<-'INH'
tab$annoType[tab$cluster %in% c(34)]<-'HTR3A'

tab$annotation[tab$cluster %in% c(26)]<-'INH'
tab$annoType[tab$cluster %in% c(26)]<-'VIP'

tab$annotation[tab$cluster %in% c(44)]<-'INH'
tab$annoType[tab$cluster %in% c(44)]<-'CORT'

######NNC
tab$annotation[tab$cluster %in% c(51,30)]<-'NNC'
tab$annoType[tab$cluster %in% c(51,30)]<-'CP'

tab$annotation[tab$cluster %in% c(1,7)]<-'NNC'
tab$annoType[tab$cluster %in% c(1,7)]<-'Oligo'

tab$annotation[tab$cluster %in% c(31)]<-'NNC'
tab$annoType[tab$cluster %in% c(31)]<-'Oligo'

tab$annotation[tab$cluster %in% c(50)]<-'NNC'
tab$annoType[tab$cluster %in% c(50)]<-'OPC'

tab$annotation[tab$cluster %in% c(22)]<-'NNC'
tab$annoType[tab$cluster %in% c(22)]<-'OPC'

tab$annotation[tab$cluster %in% c(6,17)]<-'NNC'
tab$annoType[tab$cluster %in% c(6,17)]<-'Astro'

tab$annotation[tab$cluster %in% c(37)]<-'NNC'
tab$annoType[tab$cluster %in% c(37)]<-'Ependy'

tab$annotation[tab$cluster %in% c(2)]<-'NNC'
tab$annoType[tab$cluster %in% c(2)]<-'Micro'

tab$annotation[tab$cluster %in% c(24)]<-'NNC'
tab$annoType[tab$cluster %in% c(24)]<-'Tcell.Macro'

tab$annotation[tab$cluster %in% c(27)]<-'NNC'
tab$annoType[tab$cluster %in% c(27)]<-'Vascular'

tab$annotation[tab$cluster %in% c(16)]<-'NNC'
tab$annoType[tab$cluster %in% c(16)]<-'Vascular'

tab$annotation[tab$cluster %in% c(35)]<-'NNC'
tab$annoType[tab$cluster %in% c(35)]<-'Vascular'

sce$cellType_collapsed<-tab$annoType[match(sce$k_5_louvain_1.2,tab$cluster)]


library(SingleCellExperiment)
library(ggplot2)
library(dplyr)

# Convert SingleCellExperiment to dataframe
df <- as.data.frame(colData(sce))

# Calculate proportions
df <- df %>%
    group_by(cluster_harmony, brnum) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup() %>%
    group_by(cluster_harmony) %>%
    mutate(total = sum(proportion)) %>%
    mutate(proportion = proportion / total)

# Create barplot
ggplot(df, aes(x = cluster_harmony, y = proportion, fill = as.factor(brnum))) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    ylab("Proportion") +
    xlab("Cell Type") +
    labs(fill = "brnum")





paramSweep_v3(seu, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

# First, calculate column sums
column_sums <- colSums(heat_cor)

# Then, select columns where sum is greater than 0
selected_columns <- heat_cor[, column_sums > 0]
# Define function
max_rowname_of_each_column <- function(df) {
    # Make sure df is a data frame
    if (!is.data.frame(df)) {
        stop("Input must be a dataframe")
    }

    # Initialize a vector to store the row names
    max_rownames <- character(ncol(df))

    # Loop through each column and find max
    for (i in 1:ncol(df)) {
        max_rownames[i] <- rownames(df)[which.max(df[[i]])]
    }

    # Return max rownames
    return(max_rownames)
}


features=rownames(marks2[[3]][1:25,])
plotDots(spe,group='cluster',features=c('SYT1','CHN1','MBP','GFAP','VIM','CLDN5'),color=c('white','blue')) +
    #scale_y_discrete(limits=rev(features)) +
    #scale_x_discrete(limits=rev(levels(spe$cluster)))+
    theme(axis.text.x = element_text(angle = 90,vjust=0.75),
          text = element_text(size = 14))+
    labs(x= "Cell type", y = "Gene")


plotGroupedHeatmap(spe,group='cluster',
                   features=c('SNAP25','SYT1','CHN1','PLP1','GFAP','VIM','CLDN5'),
                   color=magma(10),scale=T,center=T,cluster_cols=F,cluster_rows=F,zlim=c(-3,3))# +
    #scale_y_discrete(limits=rev(features))


plotGroupedHeatmapMod<-function (object, features, group, block = NULL, columns = NULL,
                               exprs_values = "logcounts", center = FALSE, scale = FALSE,
                               zlim = NULL, colour = color, color = NULL,
                               assay.type = exprs_values, annotation_col = NULL, ...)
{
    # Assume 'object' is a SingleCellExperiment object
    # Check if columns are null, if yes, assign sequential numbers
    if (is.null(colnames(object))) {
        colnames(object) <- seq_len(ncol(object))
    }

    # Subset the object to include only the 'features' of interest
    object <- object[features, ]

    # Subset again if 'columns' are specified
    if (!is.null(columns)) {
        object <- object[, columns]
    }

    # Compute group means
    group_means <- aggregateAcrossCells(object, ids=colData(object)[[group]],
                                        exprs_values = assay.type,
                                        fun="mean")

    # Create the heatmap
    pheatmap::pheatmap(group_means,
                       color = color,
                       show_rownames = TRUE,
                       show_colnames = TRUE,
                       annotation_col = annotation_col, ...)
}

plotGroupedHeatmapMod(spe,group='cluster',features=c('SYT1','CHN1','MBP','GFAP','VIM','CLDN5'),color=magma(25),show_colnames=T)


plotGroupedHeatmapMod <- function(spe, group, features, assay.type = "logcounts", color = pheatmap::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(25), show_colnames = F) {
    # Check inputs
    if (!group %in% colnames(colData(spe))) {
        stop(paste("Group", group, "not found in colData of the SingleCellExperiment object."))
    }

    if (!all(features %in% rownames(spe))) {
        stop(paste("Some or all features not found in the rownames of the SingleCellExperiment object."))
    }

    if (!assay.type %in% names(assays(spe))) {
        stop(paste("Assay type", assay.type, "not found in the assays of the SingleCellExperiment object."))
    }

    # Subset spe to features of interest
    spe <- spe[features, ]

    # Aggregate assay values (e.g. logcounts) across cells (by group)
    avg_logcounts <- aggregateAcrossCells(object = spe, ids = colData(spe)[[group]], exprs_values = assay.type, fun = "mean")

    # Convert to matrix
    avg_logcounts <- as.matrix(assays(avg_logcounts)[[1]])

    # Remove potential NA values
    avg_logcounts[is.na(avg_logcounts)] <- 0

    # Prepare annotation data
    annotation_data <- data.frame(Cluster = colnames(avg_logcounts))
    rownames(annotation_data) <- annotation_data$Cluster

    # Plot
    pheatmap::pheatmap(avg_logcounts, color = color, show_colnames = show_colnames, annotation_col = annotation_data)
}

