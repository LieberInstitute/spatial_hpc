###spotplots for amy/thal supplemental figure
library(SpatialExperiment)
library(ggspavis)
library(RColorBrewer)

###get rewritten plotVisium()) script
source(file=here::here('code','NMF','plotVisium_rewrite.R'))

##load palettes
load(file=here::here('plots','spatial_palette_final.rda'))

###load data
load(file=here::here('processed-data','06_clustering',
     'PRECAST','spe_precast_HE_domain.rda'))

spatial.palette3<-c("#BEDDBA", "#eae8e4", "#E8BBC6","#A1BAD8")
names(spatial.palette3)<-levels(speb$broad.domain)

brains<-unique(spe$brnum)

speb <- spe[, (colData(spe)$brnum == brains[[4]])]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)


pdf(file=here::here('plots','figures','supp_figures','figure_S24','slc_spotplot.pdf'),height=7,w=7)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'SLC17A6',
    highlight = "broad.domain",
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    values=spatial.palette3,
    image=F)+
scale_fill_distiller(
    type = "seq",
    palette = rev('Greys'),
    direction=1)
dev.off()

pdf(file=here::here('plots','figures','supp_figures','figure_S24','oprm1_spotplot.pdf'),height=7,w=7)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'OPRM1',
    highlight = "broad.domain",
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    values=spatial.palette3,
    image=F)+
scale_fill_distiller(
    type = "seq",
    palette = rev('Greys'),
    direction=1)
dev.off()

pdf(file=here::here('plots','figures','supp_figures','figure_S24','CACNG4_spotplot.pdf'),height=7,w=7)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'CACNG4',
    highlight = "broad.domain",
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    values=spatial.palette3,
    image=F)+
scale_fill_distiller(
    type = "seq",
    palette = rev('Greys'),
    direction=1)
dev.off()


pdf(file=here::here('plots','figures','supp_figures','figure_S24','CDH22_spotplot.pdf'),height=7,w=7)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'CDH22',
    highlight = "broad.domain",
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    values=spatial.palette3,
    image=F)+
scale_fill_distiller(
    type = "seq",
    palette = rev('Greys'),
    direction=1)
dev.off()


#####thal spotplots
speb <- spe[, (colData(spe)$brnum == brains[[7]])]
speb<-speb[,!speb$sample_id=='V11L05-335_D1']
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)


pdf(file=here::here('plots','figures','supp_figures','figure_S24','TCF7L2_spotplot.pdf'),height=7,w=7)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'TCF7L2',
    highlight = "broad.domain",
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    values=spatial.palette3,
    image=F)+
scale_fill_distiller(
    type = "seq",
    palette = rev('Greys'),
    direction=1)
dev.off()

pdf(file=here::here('plots','figures','supp_figures','figure_S24','SHOX2_spotplot.pdf'),height=7,w=7)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'SHOX2',
    highlight = "broad.domain",
    facets = "sample_id",
    assay = "logcounts",
    trans = "identity",
    x_coord = NULL,
    y_coord = NULL,
    y_reverse = TRUE,
    sample_ids = NULL,
    image_ids = NULL,
    values=spatial.palette3,
    image=F)+
scale_fill_distiller(
    type = "seq",
    palette = rev('Greys'),
    direction=1)
dev.off()
