setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
library(ggspavis)
library(SpatialExperiment)
library(here)
library(scater)
library(scran)
library(ggrastr)
load(file=here::here('processed-data','06_clustering','PRECAST','spe_precast_HE.rda'))

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
tab$annotation[tab$cluster %in% c(15)]<-'SLM.WM'
tab$annotation[tab$cluster %in% c(16)]<-'WM.2'
tab$annotation[tab$cluster %in% c(17)]<-'CA2.4.2'
tab$annotation[tab$cluster %in% c(18)]<-'WM.3'


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
tab$annotation2[tab$cluster %in% c(15)]<-'Neuropil'
tab$annotation2[tab$cluster %in% c(16)]<-'WM'
tab$annotation2[tab$cluster %in% c(17)]<-'Neuron'
tab$annotation2[tab$cluster %in% c(18)]<-'WM'

spe$cluster<-tab$annotation[match(spe$PRECAST_k18,tab$cluster)]
spe$broad<-factor(tab$annotation2[match(spe$PRECAST_k18,tab$cluster)])

spe$cluster<-tab$annotation[match(spe$cluster,tab$annotation)]
spe$cluster<-factor(spe$cluster,levels=c("GCL","CA2.4.1", "CA2.4.2",
                                                       "CA1.1", "CA1.2","SUB", "SUB.RHP",
                                                       "RHP","GABA","SL.SR","ML", "SR.SLM","SLM.WM",
                                                       paste0("WM.",c(1,2,3)),"Vascular","Choroid"))

spe$domain<-spe$cluster
levels(spe$domain)[2]<-'CA2.4'
levels(spe$domain)[3]<-'CA2.4'
levels(spe$domain)[3]<-'CA1'
levels(spe$domain)[4]<-'CA1'

clustering_name <- 'domain'




#make sure these are factors
spe$brnum<-factor(spe$brnum)
spe$sample_id<-factor(spe$sample_id)

#make clustering plots
brains <- unique(spe$brnum)
k<-16
p<-list()
for (i in seq_along(clustering_name)) {
    p[[i]]<-list()
    for (j in seq_along(brains)) {
        speb <- spe[, (colData(spe)$brnum == brains[j])]
        speb$sample_id <- droplevels(speb$sample_id)
        speb$sample_id <- as.character(speb$sample_id)
        samples <- unique(speb$sample_id)
        speb$sample_id <- factor(speb$sample_id, levels = samples)
        samples
        speb$brnum <- droplevels(speb$brnum)

   palette=spatial.palette



        names(palette) <- levels(speb[[clustering_name[i]]])
        print(paste0("printing plot",i,"_",j))
        p[[i]][[j]]<-plotVisium(
            speb,
            spots = TRUE,
            fill = clustering_name[i],
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
            palette = palette
        )+ ggplot2::theme(legend.text = ggplot2::element_text(size = 12),
                           plot.title = ggplot2::element_text(size = 18))
    }}

for (i in seq_along(clustering_name)) {
    pdf(
        file = here::here("plots", "06_clustering", "PRECAST",
                          paste0(clustering_name[i], ".pdf")))


    for (j in seq_along(p[[i]])) {
        print(p[[i]][[j]])
    }

    dev.off()
}


#make marker gene plots
speb <- spe[, (colData(spe)$sample_id == 'V11L05-333_D1')]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

pdf(here::here('plots','figures','figure_2','visium_ppfia2.pdf'),h=4,w=3)
plotVisium(
  speb,
  spots = TRUE,
  fill = 'PPFIA2',
  highlight = NULL,
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
  highlight = NULL,
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
for (i in seq_along(clustering_name)) {
    p[[i]]<-list()
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

        print(paste0("printing plot",i,"_",j))
        p[[i]][[j]]<-plotVisium(
            speb,
            spots = TRUE,
            fill = clustering_name[i],
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
    }}

pdf(
        file = here::here("plots", "figures", "supp_figures",
                          paste0(clustering_name[i] ,".pdf")))


    for (j in seq_along(p[[i]])) {
        p[[i]][[j]]<-rasterize(p[[i]][[j]],dpi=350,layers='Point')
        print(p[[i]][[j]])
    }

dev.off()


#####nnSVG plots######
speb <- spe[, (colData(spe)$brnum == brains[2])]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

pdf(
    file = here::here("plots", "figures", "supp_figures",
                      "nnSVG_top_genes.pdf"))
genes<-c('MBP','GFAP','MT_CO3','FOLR1','UCHL1','KRT17')
for(i in 1:6){
p<-rasterize(
    plotVisium(
        speb,
        spots = TRUE,
        fill = genes[i],
        highlight = NULL,
        facets = "sample_id",
        assay = "logcounts",
        trans = "identity",
        x_coord = NULL,
        y_coord = NULL,
        y_reverse = TRUE,
        sample_ids = NULL,
        image_ids = NULL,
        #palette = spatial.palette,
        image=F)+
        scale_fill_distiller(
            type = "seq",
            palette = rev('Greys'),
            direction=1)+theme(legend.text=element_text(size = 12),
                               legend.title=element_text(size = 20))
,dpi=200,layers='Point')
print(p)
}
dev.off()
