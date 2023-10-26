#####supplementary figs related to figs 1 and 2
setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("ggspavis"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("ggrastr"))
suppressPackageStartupMessages(library("scater"))



load(here("processed-data", "04_QC", "spe_QC_allSamples.Rdata"), verbose = TRUE)

##Histology fig
pdf(here::here('plots','figures','supp_figures','visium_histology.pdf'),h=4,w=4)
brains <- unique(spe$brnum)

for(i in 1:length(brains)){
print(paste0('printing plot ',i))
speb <- spe[, (colData(spe)$brnum == brains[[i]])]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)
print(plotVisium(
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
    image=T))
}
dev.off()
rm(spe,speb)

load(file=here::here('processed-data','05_preprocess_batchCorrection','spe_norm.rda'))
#SNAP25
pdf(here::here('plots','figures','supp_figures','visium_SNAP25.pdf'),h=7,w=7)
brains <- unique(spe$brnum)

for(i in 1:length(brains)){
    print(paste0('printing plot ',i))
    speb <- spe[, (colData(spe)$brnum == brains[[i]])]
    speb$sample_id <- droplevels(speb$sample_id)
    speb$sample_id <- as.character(speb$sample_id)
    samples <- unique(speb$sample_id)
    speb$sample_id <- factor(speb$sample_id, levels = samples)
    samples
    speb$brnum <- droplevels(speb$brnum)
    print(
        rasterize(
        plotVisium(
        speb,
        spots = T,
        fill = 'SNAP25',
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
        image=F)+scale_fill_distiller(
            type = "seq",
            palette = rev('Greys'),
            direction=1),dpi=250))
}
dev.off()

#MBP
pdf(here::here('plots','figures','supp_figures','visium_MBP.pdf'),h=7,w=7)
brains <- unique(spe$brnum)

for(i in 1:length(brains)){
    print(paste0('printing plot ',i))
    speb <- spe[, (colData(spe)$brnum == brains[[i]])]
    speb$sample_id <- droplevels(speb$sample_id)
    speb$sample_id <- as.character(speb$sample_id)
    samples <- unique(speb$sample_id)
    speb$sample_id <- factor(speb$sample_id, levels = samples)
    samples
    speb$brnum <- droplevels(speb$brnum)
    print(
        rasterize(
            plotVisium(
        speb,
        spots = T,
        fill = 'MBP',
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
        image=F)+scale_fill_distiller(
            type = "seq",
            palette = rev('Greys'),
            direction=1),dpi=250))
}
dev.off()

