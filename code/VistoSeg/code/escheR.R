setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
  library("here")
  library("SpatialExperiment")
  library("spatialLIBD")
  library("rtracklayer")
  library("lobstr")
  library("sessioninfo")
  library("escheR")
})
load('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/02_build_spe/spe_bayes_precast.rda')

sample_id = c("V11A20-297_A1", "V11A20-297_B1", "V11A20-297_C1", "V11A20-297_D1", "V11L05-333_A1", "V11L05-333_B1", "V11L05-333_C1", "V11L05-333_D1") 
spea <- spe[, spe$sample_id %in% sample_id]

segmentations_list <-
  lapply(sample_id, function(sampleid) {
    file <-
      here(
        "processed-data",
        "01_spaceranger",
        "spaceranger_2022-04-12_SPag033122",
        sampleid,
        "outs",
        "spatial",
        "tissue_spot_counts.csv"
      )
    if (!file.exists(file)) {
      return(NULL)
    }
    x <- read.csv(file)
    x$key <- paste0(x$barcode, "_", sampleid)
    return(x)
  })

## Merge them (once the these files are done, this could be replaced by an rbind)
segmentations <-
  Reduce(function(...) {
    merge(..., all = TRUE)
  }, segmentations_list[lengths(segmentations_list) > 0])

## Add the information
segmentation_match <- match(spea$key, segmentations$key)
segmentation_info <-
  segmentations[segmentation_match, -which(
    colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
  )]
colData(spea) <- cbind(colData(spea), segmentation_info)
colData(spea)$sample_id = as.character(colData(spea)$sample_id)

pdf(here("plots", "escheR.pdf"), width = 10, height = 10)
speb = spea[, which(spea$sample_id == sample_id[1])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")

speb = spea[, which(spea$sample_id == sample_id[2])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")

speb = spea[, which(spea$sample_id == sample_id[3])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")

speb = spea[, which(spea$sample_id == sample_id[4])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")

speb = spea[, which(spea$sample_id == sample_id[5])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")

speb = spea[, which(spea$sample_id == sample_id[6])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")

speb = spea[, which(spea$sample_id == sample_id[7])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")

speb = spea[, which(spea$sample_id == sample_id[8])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")

dev.off()

##
pdf(here("plots", "escheR_swap.pdf"), width = 10, height = 10)
speb = spea[, which(spea$sample_id == sample_id[1])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "Pmask_dark_blue")
p1 |> add_fill(var = "PRECAST_k16_nnSVG")+scale_fill_viridis_d(option = "H")

speb = spea[, which(spea$sample_id == sample_id[2])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "Pmask_dark_blue")
p1 |> add_fill(var = "PRECAST_k16_nnSVG")

speb = spea[, which(spea$sample_id == sample_id[3])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "Pmask_dark_blue")
p1 |> add_fill(var = "PRECAST_k16_nnSVG")

speb = spea[, which(spea$sample_id == sample_id[4])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "Pmask_dark_blue")
p1 |> add_fill(var = "Pmask_dark_blue")

speb = spea[, which(spea$sample_id == sample_id[5])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "Pmask_dark_blue")
p1 |> add_fill(var = "PRECAST_k16_nnSVG")

speb = spea[, which(spea$sample_id == sample_id[6])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "Pmask_dark_blue")
p1 |> add_fill(var = "PRECAST_k16_nnSVG")

speb = spea[, which(spea$sample_id == sample_id[7])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "Pmask_dark_blue")
p1 |> add_fill(var = "PRECAST_k16_nnSVG")

speb = spea[, which(spea$sample_id == sample_id[8])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "Pmask_dark_blue")
p1 |> add_fill(var = "PRECAST_k16_nnSVG")

dev.off()



## 
pdf(here("plots", "nuclei density plots.pdf"), width = 10, height = 10)
for (i in seq_along(sample_id)){
 speb = spea[, which(spea$sample_id == sample_id[i])]
 den <- density(speb$Pmask_dark_blue)
 plot(den, frame = FALSE, col = "blue",main = sample_id[i])
}
dev.off()

colData(spea)$nuc_bin = 0
colData(spea)$nuc_bin[colData(spea)$Pmask_dark_blue>0.05] = 1
colData(spea)$nuc_bin[colData(spea)$Pmask_dark_blue>0.25] = 2
colData(spea)$nuc_bin = as.factor(colData(spea)$nuc_bin)
spea$nuc_bin = factor(spea$nuc_bin, levels = c("1" ,"2"))

pdf(here("plots", "escheR_3bin_emptyspots.pdf"), width = 10, height = 10)
speb = spea[, which(spea$sample_id == sample_id[1])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")+scale_colour_viridis_d(option = "H")
p1 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 16))

speb = spea[, which(spea$sample_id == sample_id[2])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")+scale_colour_viridis_d(option = "H")
p1 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 16))

speb = spea[, which(spea$sample_id == sample_id[3])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")+scale_colour_viridis_d(option = "H")
p1 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 16))

speb = spea[, which(spea$sample_id == sample_id[4])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")+scale_colour_viridis_d(option = "H")
p1 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 16))

speb = spea[, which(spea$sample_id == sample_id[5])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")+scale_colour_viridis_d(option = "H")
p1 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 16))

speb = spea[, which(spea$sample_id == sample_id[6])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")+scale_colour_viridis_d(option = "H")
p1 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 16))

speb = spea[, which(spea$sample_id == sample_id[7])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")+scale_colour_viridis_d(option = "H")
p1 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 16))

speb = spea[, which(spea$sample_id == sample_id[8])]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")+scale_colour_viridis_d(option = "H")
p1 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 16))

dev.off()
