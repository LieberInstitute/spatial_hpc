#####plots for thal/amy ID in snRNAseq data
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(here)

sce_neuron<-sce[,sce$broad.cell.class=='Neuron']
#####thal marker genes
pdf(file=here::here('plots','figures','supp_figures','figure_S24','tcf7l2.pdf'),height=2,w=4.2)
ggcells(sce_neuron,
        mapping=aes(x=superfine.cell.class, y=TCF7L2,fill=fine.cell.class)) +
        geom_boxplot(outlier.size=0.1)+
        theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 10,colour='black'),
                         legend.position = 'bottom',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))+
        scale_fill_manual(values=sn.fine.palette)
dev.off()

pdf(file=here::here('plots','figures','supp_figures','figure_S24','shox2.pdf'),height=2,w=4.2)
ggcells(sce_neuron,
        mapping=aes(x=superfine.cell.class, y=SHOX2,fill=fine.cell.class)) +
        geom_boxplot(outlier.size=0.1)+
        theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 10,colour='black'),
                         legend.position = 'none',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))+
        scale_fill_manual(values=sn.fine.palette)
dev.off()

pdf(file=here::here('plots','figures','supp_figures','figure_S24','slc17a6.pdf'),height=2,w=4.2)
ggcells(sce_neuron,
        mapping=aes(x=superfine.cell.class, y=SLC17A6,fill=fine.cell.class)) +
        geom_boxplot(outlier.size=0.1)+
        theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 10,colour='black'),
                         legend.position = 'none',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))+
        scale_fill_manual(values=sn.fine.palette)
dev.off()

pdf(file=here::here('plots','figures','supp_figures','figure_S24','ptger3.pdf'),height=2,w=4.2)
ggcells(sce_neuron,
        mapping=aes(x=superfine.cell.class, y=PTGER3,fill=fine.cell.class)) +
        geom_boxplot(outlier.size=0.1)+
        theme(axis.text.x = element_text(angle = 90),
                         text=element_text(size = 10,colour='black'),
                         legend.position = 'none',
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.line = element_line(colour = "black"))+
        scale_fill_manual(values=sn.fine.palette)
dev.off()

###spotplots
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

pdf(file=here::here('plots','figures','supp_figures','figure_S24','SCN11A_spotplot.pdf'),height=7,w=7)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'SCN11A',
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

pdf(file=here::here('plots','figures','supp_figures','figure_S24','DCSTAMP_spotplot.pdf'),height=7,w=7)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'DCSTAMP',
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

pdf(file=here::here('plots','figures','supp_figures','figure_S24','M1AP_spotplot.pdf'),height=7,w=7)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'M1AP',
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

pdf(file=here::here('plots','figures','supp_figures','figure_S24','PAPPA2_spotplot.pdf'),height=7,w=7)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'PAPPA2',
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

pdf(file=here::here('plots','figures','supp_figures','figure_S24','PENK_spotplot.pdf'),height=7,w=7)
plotVisium(
    speb,
    spots = TRUE,
    fill = 'PENK',
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

















