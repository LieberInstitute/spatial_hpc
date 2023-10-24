setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("ggspavis"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("ggrastr"))

load(here("processed-data", "04_QC", "spe_QC_allSamples.Rdata"), verbose = TRUE)
spe$discard_auto_br <- spe$low_sum_br | spe$low_detected_br
spe$discard_auto_id <- spe$low_sum_id | spe$low_detected_id

# ## low sum/detected genes
pdf(file=here::here('plots','figures','supp_figures','QC_plots_sampleID.pdf'))
p1<-plotColData(spe, x = "sample_id", y = "sum", colour_by = "low_sum_id",point_alpha=1,point_size=0.25) +
  scale_y_log10() +
  ggtitle("Total UMIs") +theme(axis.text.x = element_text(size = 10, angle = 90))
p2<-plotColData(spe, x = "sample_id", y = "detected", colour_by = "low_detected_id",point_alpha=1,point_size=0.25) +
  scale_y_log10() +
  ggtitle("Detected genes") +theme(axis.text.x = element_text(size = 10, angle = 90))
grid.arrange(p1, p2,nrow=2)
dev.off()
#  facet_wrap(~ sce$round, scales = "free_x", nrow = 1)
# +
#   geom_hline(yintercept = 1000) ## hline doesn't work w/ facet_wrap?

pdf(here::here('plots','figures','supp_figures','mito_rate.pdf'),h=7,w=7)
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
    fill = 'mito_percent',
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
      direction=1,
      limits=c(0,50)),dpi=500))
}
dev.off()

pdf(here::here('plots','figures','supp_figures','discarded_spots.pdf'),h=7,w=7)
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
    fill = 'discard_auto_id',
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
    image=F)+scale_fill_manual(values=c('#dddddd','magenta'))),dpi=500)#+scale_fill_distiller(
  #type = "seq",
  # palette = rev('RdPu'),
  # direction=1))
}
dev.off()
