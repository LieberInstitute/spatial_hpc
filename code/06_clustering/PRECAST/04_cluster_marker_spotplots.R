setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
library(ggspavis)
library(SpatialExperiment)
library(here)
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

clustering_name <- 'cluster'

#make sure these are factors
spe$brnum<-factor(spe$brnum)
spe$sample_id<-factor(spe$sample_id)

#make clustering plots
brains <- unique(spe$brnum)
k<-18
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

        palette <- c(
          "#008000",  # Level 1 Brighter Forest Green
          "#FF1493",  # Level 2 Bright Hot Pink
          "#FFB6C1",  # Level 3 Medium Pink
          "#800080",  # Level 4 Purple (darker muted purple)
          "#D8BFD8",  # Level 5 Thistle (light muted purple)
          "#FFEC8B",  # Level 6 Light Goldenrod
          "#DAA520",  # Level 7 Goldenrod
          "#FF0000",  # Level 8 Bright Red (less brown, more red)
          "#98FF98",  # Level 9 Mint Green (mintier)
          "#111111",  # Level 10 Lighter Black
          "#dddddd",  # Level 11 Lighter Gray
          "#A9A9A9",  # Level 12 Medium Gray (lighter)
          "#000080",  # Level 13 Lighter Deep Navy
          "#87CEEB",  # Level 14 Sky Blue
          "#4682B4",  # Level 15 Steel Blue
          "#3333DD",  # Level 16 Lighter Dark Blue
          "#F4A460",  # Level 17 Sandy Brown
          "#8B4513"   # Level 18 Saddle Brown
        )
        
        
    
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
