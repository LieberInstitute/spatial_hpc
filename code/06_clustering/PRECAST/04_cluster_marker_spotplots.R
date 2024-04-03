
setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
library(ggspavis)
library(SpatialExperiment)
library(here)
library(scater)
library(scran)
library(ggrastr)

###load data
load(file=here::here('processed-data','06_clustering','PRECAST','spe_precast_HE.rda'))

###get rewritten plotVisium()) script
source(file=here::here('code','NMF','plotVisium_rewrite.R'))

##load palettes
load(file=here::here('plots','spatial_palette_final.rda'))

####make annotation data frame for cluster
tab<-data.frame('cluster'=1:18,'annotation2'=rep(NA,18),'annotation'=rep(NA,18))
tab$annotation[tab$cluster %in% c(4)]<-'GCL'
tab$annotation[tab$cluster %in% c(1)]<-'SUB.RHP'
tab$annotation[tab$cluster %in% c(2)]<-'ML'
tab$annotation[tab$cluster %in% c(3)]<-'CA1.1'
tab$annotation[tab$cluster %in% c(5)]<-'SR.SLM'
tab$annotation[tab$cluster %in% c(6)]<-'WM.1'
tab$annotation[tab$cluster %in% c(7)]<-'CA2.4.1'
tab$annotation[tab$cluster %in% c(8)]<-'RHP'
tab$annotation[tab$cluster %in% c(9)]<-'CA1.2'
tab$annotation[tab$cluster %in% c(10)]<-'GABA'
tab$annotation[tab$cluster %in% c(11)]<-'Vascular'
tab$annotation[tab$cluster %in% c(12)]<-'Choroid'
tab$annotation[tab$cluster %in% c(13)]<-'SL.SR'
tab$annotation[tab$cluster %in% c(14)]<-'SUB'
tab$annotation[tab$cluster %in% c(15)]<-'SLM.SGZ'
tab$annotation[tab$cluster %in% c(16)]<-'WM.2'
tab$annotation[tab$cluster %in% c(17)]<-'CA2.4.2'
tab$annotation[tab$cluster %in% c(18)]<-'WM.3'


###make annotation data frame for broad domains
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
tab$annotation2[tab$cluster %in% c(11)]<-'Vasc_CSF'
tab$annotation2[tab$cluster %in% c(12)]<-'Vasc_CSF'
tab$annotation2[tab$cluster %in% c(13)]<-'Neuropil'
tab$annotation2[tab$cluster %in% c(14)]<-'Neuron'
tab$annotation2[tab$cluster %in% c(15)]<-'Neuropil'
tab$annotation2[tab$cluster %in% c(16)]<-'WM'
tab$annotation2[tab$cluster %in% c(17)]<-'Neuron'
tab$annotation2[tab$cluster %in% c(18)]<-'WM'

###make cluster and domain columns in spe
spe$cluster<-tab$annotation[match(spe$PRECAST_k18,tab$cluster)]
spe$broad.domain<-factor(tab$annotation2[match(spe$PRECAST_k18,tab$cluster)],levels=c('Neuron','Neuropil','WM','Vasc_CSF'))


###make cluster a factor
spe$cluster<-factor(spe$cluster,levels=c("GCL","CA2.4.1", "CA2.4.2",
                                                       "CA1.1", "CA1.2","SUB", "SUB.RHP",
                                                       "RHP","GABA","SL.SR","ML", "SR.SLM","SLM.SGZ",
                                                       paste0("WM.",c(1,2,3)),"Vascular","Choroid"))

###group CA2-4 and CA1 clusters
spe$domain<-spe$cluster
levels(spe$domain)[2]<-'CA2.4'
levels(spe$domain)[3]<-'CA2.4'
levels(spe$domain)[3]<-'CA1'
levels(spe$domain)[4]<-'CA1'

####save spe with labeled domains
save(spe,file=here::here('processed-data','06_clustering','PRECAST','spe_precast_HE_domain.rda'))


############make spotplots#####################
clustering_name <- 'domain'

#make sure these are factors
spe$brnum<-factor(spe$brnum)
spe$sample_id<-factor(spe$sample_id)

#make clustering plots
brains <- unique(spe$brnum)
k<-18
p<-list()
    for (j in seq_along(brains)) {
        speb <- spe[, (colData(spe)$brnum == brains[j])]
        speb$sample_id <- droplevels(speb$sample_id)
        speb$sample_id <- as.character(speb$sample_id)
        samples <- unique(speb$sample_id)
        speb$sample_id <- factor(speb$sample_id, levels = samples)
        samples
        speb$brnum <- droplevels(speb$brnum)

   palette=srt.palette



        names(palette) <- levels(speb[[clustering_name]])
        print(paste0("printing plot",j))
        p[[j]]<-plotVisium(
            speb,
            spots = TRUE,
	    highlight='domain',
	    values=palette,
            fill = clustering_name,
            #highlight = NULL,
            facets = "sample_id",
            image = FALSE,
            assay = "logcounts",
            trans = "identity",
            x_coord = NULL,
            y_coord = NULL,
            y_reverse = TRUE,
            sample_ids = NULL,
            image_ids = NULL,
            palette = palette
        )+ ggplot2::theme(legend.text = ggplot2::element_text(size = 12),
                           plot.title = ggplot2::element_text(size = 18))
    }


pdf(file = here::here("plots", "06_clustering", "PRECAST",
               paste0(clustering_name, ".pdf")))

    for (j in seq_along(p)) {
        print(p[[j]]) }

dev.off()



############make marker gene plots for fig 2############
speb <- spe[, (colData(spe)$sample_id == 'V11L05-333_D1')]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

###make highlight palette
spatial.palette3<-c("#BEDDBA", "#eae8e4", "#E8BBC6","#A1BAD8")
names(spatial.palette3)<-levels(speb$broad.domain)
pdf(here::here('plots','figures','figure_2','visium_ppfia2.pdf'),h=4,w=3)
plotVisium(
  speb,
  spots = TRUE,
  fill = 'PPFIA2',
  highlight = 'broad.domain',
  values=spatial.palette3,
  facets = NULL,
  assay = "logcounts",
  trans = "identity",
  x_coord = NULL,
  y_coord = NULL,
  y_reverse = TRUE,
  sample_ids = NULL,
  image_ids = NULL,
  image=F)+scale_fill_distiller(
    type = "seq",
    palette = rev('Greys'),
    direction=1)
dev.off()

pdf(here::here('plots','figures','figure_2','visium_prkcg.pdf'),h=4,w=3)
plotVisium(
  speb,
  spots = TRUE,
  fill = 'PRKCG',
  highlight = 'broad.domain',
  values=spatial.palette3,
  facets = NULL,
  assay = "logcounts",
  trans = "identity",
  x_coord = NULL,
  y_coord = NULL,
  y_reverse = TRUE,
  sample_ids = NULL,
  image_ids = NULL,
  #palette = as.vector(palette36.colors(18)),
  image=F
)+scale_fill_distiller(
  type = "seq",
  palette = rev('Greys'),
  direction=1)
dev.off()

pdf(here::here('plots','figures','figure_2','visium_APOC1.pdf'),h=4,w=3)
plotVisium(
  speb,
  spots = TRUE,
  fill = 'APOC1',
  highlight = 'broad.domain',
  values=spatial.palette3,
  facets = "sample_id",
  assay = "logcounts",
  trans = "identity",
  x_coord = NULL,
  y_coord = NULL,
  y_reverse = TRUE,
  sample_ids = NULL,
  image_ids = NULL,
  #palette = as.vector(palette36.colors(18)),
  image=F
)+scale_fill_distiller(
  type = "seq",
  palette = rev('Greys'),
  direction=1)
dev.off()

pdf(here::here('plots','figures','figure_2','visium_SFRP2.pdf'),h=4,w=3)
plotVisium(
  speb,
  spots = TRUE,
  fill = 'SFRP2',
  highlight = 'broad.domain',
  values=spatial.palette3,
  facets = "sample_id",
  assay = "logcounts",
  trans = "identity",
  x_coord = NULL,
  y_coord = NULL,
  y_reverse = TRUE,
  sample_ids = NULL,
  image_ids = NULL,
  #palette = as.vector(palette36.colors(18)),
  image=F
)+scale_fill_distiller(
  type = "seq",
  palette = rev('Greys'),
  direction=1)
dev.off()


###make QC metric plots for fig s17:
pdf(file=here::here('plots','06_clustering','PRECAST','violin_plots_QC_metrics.pdf'),h=8.5,w=8.5)
p1<-plotColData(spe,x='domain',y='sum',point_size=0.5,add_legend=F,colour_by='domain')+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
               width = 0.3)+
  theme(axis.text.x = element_text(angle = 90),
        text=element_text(size = 10))+
        ggtitle('Library size')+
        labs(y='library size')+
        scale_y_log10()+scale_color_manual(values=spatial.palette)
p1<-rasterize(p1,dpi = 400,layers='Point')

p2<-plotColData(spe,x='domain',y='detected',point_size=0.5,add_legend=F,colour_by='domain')+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
               width = 0.3)+
  theme(axis.text.x = element_text(angle = 90),
        text=element_text(size = 10))+
  ggtitle('Detected genes')+
  labs(y='detected genes')+
  scale_y_log10()+scale_color_manual(values=spatial.palette)
p2<-rasterize(p2,dpi = 350,layers='Point')

p3<-plotColData(spe,x='domain',y='subsets_Mito_percent',point_size=0.5,add_legend=F,colour_by='domain')+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
               width = 0.3)+
  theme(axis.text.x = element_text(angle = 90),
        text=element_text(size = 10))+
  ggtitle('Mitochondrial gene expression percentage')+
  labs(y='mito expression percent')+
  scale_color_manual(values=spatial.palette)
p3<-rasterize(p3,dpi = 350,layers='Point')

  grid.arrange(p1,p2,p3,ncol=1)

dev.off()

###Manual annotation plots
###some spots are NA lets remove them
spe<-spe[,!is.na(spe$ManualAnnotation)]
#make clustering plots
brains <- unique(spe$brnum)
k<-15
clustering_name<-'ManualAnnotation'
p<-list()
    for (j in seq_along(brains)) {
        speb <- spe[, (colData(spe)$brnum == brains[j])]
        speb$sample_id <- droplevels(speb$sample_id)
        speb$sample_id <- as.character(speb$sample_id)
        samples <- unique(speb$sample_id)
        speb$sample_id <- factor(speb$sample_id, levels = samples)
        samples
        speb$brnum <- droplevels(speb$brnum)
        palette<-Polychrome::palette36.colors(15)
        names(palette)<-levels(spe$ManualAnnotation)

        print(paste0("printing plot",j))
        p[[j]]<-plotVisium(
            speb,
            spots = TRUE,
            fill = clustering_name,
            highlight = NULL,
            facets = "sample_id",
            image = FALSE,
            assay = "logcounts",
            trans = "identity",
            x_coord = NULL,
            y_coord = NULL,
            y_reverse = TRUE,
            sample_ids = NULL,
            image_ids = NULL,
            palette=palette
        )+ ggplot2::theme(legend.text = ggplot2::element_text(size = 12),
                          plot.title = ggplot2::element_text(size = 18))
    }

pdf(
        file = here::here("plots", "figures", "supp_figures",
                          paste0(clustering_name[i] ,".pdf")))


    for (j in seq_along(p)) {
        p[[j]]<-rasterize(p[[j]],dpi=350,layers='Point')
        print(p[[j]])
    }

dev.off()

#####heatmaps for precast k18 grouping justification
pdf(file=here::here('plots','figures','supp_figures','figure_S12','heatmap.pdf'))
plotGroupedHeatmap(spe,group = 'cluster',features=c('FIBCD1','FNDC1','COL5A2','TSPAN18','AMPH','KCNQ5'),
scale=T,center=T,clustering_distance_cols='canberra',cluster_rows=F)
dev.off()

pdf(file=here::here('plots','figures','supp_figures','figure_S12','heatmap_astro.pdf'))
plotGroupedHeatmap(spe,group = 'cluster',features=c('GFAP','AQP4','APOE','ID4','SOX2'),
scale=T,center=T,clustering_distance_cols='canberra',cluster_rows=F)
dev.off()

###boxplots showing differences in nuclei numbers
spe<-spe[,spe$cluster %in% c('CA2.4.1','CA2.4.2','CA1.1','CA1.2')]
spe$cluster<-droplevels(spe$cluster)

pdf(file=here::here('plots','figures','supp_figures','figure_S12','boxplots_nuclei_count.pdf'),h=1.5,w=3.75)
ggcells(spe,
        mapping=aes(x=cluster, y=count,fill=cluster)) +
        geom_boxplot(outliers=FALSE)+
        theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 10,colour='black'),
                         legend.position = 'none',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))+
	labs(x='PRECAST clusters (k=18)',y='number of nuclei per spot')
dev.off()

pdf(file=here::here('plots','figures','supp_figures','figure_S12','boxplots_libsize.pdf'),h=2,w=3.75)
ggcells(spe,
        mapping=aes(x=cluster, y=sum,fill=cluster)) +
        geom_boxplot(outliers=FALSE)+
        theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 10,colour='black'),
                         legend.position = 'none',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))+
	labs(x='PRECAST clusters (k=18)',y='library size')
dev.off()
