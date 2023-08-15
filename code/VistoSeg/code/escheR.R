setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
  library("here")
  library("SpatialExperiment")
  library("spatialLIBD")
  library("rtracklayer")
  library("lobstr")
  library("sessioninfo")
  library("escheR")
  library("ggplot2")
  library("gridExtra")
})

load(here("processed-data","02_build_spe","spe_nmf_final.rda"))
t = colData(spe)[colData(spe)$slide=="V12F14-051",]
temp = paste0(sapply(strsplit(t$key,"Br"),'[',1),t$sample_id)
colData(spe)$key[spe$slide=="V12F14-051"]=temp

spaceranger_dirs = read.csv(file.path(here::here("code","VistoSeg","code","samples.txt")), header = FALSE, sep = '\t', stringsAsFactors = FALSE, col.names = c('SPpath','sample_id','brain'))
spaceranger_dirs = spaceranger_dirs[1:36,]
spaceranger_dirs$SPpath = paste0(spaceranger_dirs$SPpath,"outs/spatial/tissue_spot_counts.csv")

segmentations_list <-
  lapply(spaceranger_dirs$sample_id, function(sampleid) {
    file <-spaceranger_dirs$SPpath[spaceranger_dirs$sample_id == sampleid]
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
segmentation_match <- match(spe$key, segmentations$key)
segmentation_info <-
  segmentations[segmentation_match, -which(
    colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
  )]
colData(spe) <- cbind(colData(spe), segmentation_info)
colData(spe)$sample_id = as.character(colData(spe)$sample_id)

# pdf(here("plots", "escheR.pdf"), width = 10, height = 10)
# speb = spea[, which(spea$sample_id == sample_id[1])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
# p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")
# 
# speb = spea[, which(spea$sample_id == sample_id[2])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
# p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")
# 
# speb = spea[, which(spea$sample_id == sample_id[3])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
# p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")
# 
# speb = spea[, which(spea$sample_id == sample_id[4])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
# p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")
# 
# speb = spea[, which(spea$sample_id == sample_id[5])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
# p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")
# 
# speb = spea[, which(spea$sample_id == sample_id[6])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
# p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")
# 
# speb = spea[, which(spea$sample_id == sample_id[7])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
# p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")
# 
# speb = spea[, which(spea$sample_id == sample_id[8])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "PRECAST_k16_nnSVG")
# p1 |> add_fill(var = "Pmask_dark_blue")+scale_fill_viridis_c(option = "H")
# 
# dev.off()
# 
# ##
# pdf(here("plots", "escheR_swap.pdf"), width = 10, height = 10)
# speb = spea[, which(spea$sample_id == sample_id[1])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "Pmask_dark_blue")
# p1 |> add_fill(var = "PRECAST_k16_nnSVG")+scale_fill_viridis_d(option = "H")
# 
# speb = spea[, which(spea$sample_id == sample_id[2])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "Pmask_dark_blue")
# p1 |> add_fill(var = "PRECAST_k16_nnSVG")
# 
# speb = spea[, which(spea$sample_id == sample_id[3])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "Pmask_dark_blue")
# p1 |> add_fill(var = "PRECAST_k16_nnSVG")
# 
# speb = spea[, which(spea$sample_id == sample_id[4])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "Pmask_dark_blue")
# p1 |> add_fill(var = "Pmask_dark_blue")
# 
# speb = spea[, which(spea$sample_id == sample_id[5])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "Pmask_dark_blue")
# p1 |> add_fill(var = "PRECAST_k16_nnSVG")
# 
# speb = spea[, which(spea$sample_id == sample_id[6])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "Pmask_dark_blue")
# p1 |> add_fill(var = "PRECAST_k16_nnSVG")
# 
# speb = spea[, which(spea$sample_id == sample_id[7])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "Pmask_dark_blue")
# p1 |> add_fill(var = "PRECAST_k16_nnSVG")
# 
# speb = spea[, which(spea$sample_id == sample_id[8])]
# p <- make_escheR(speb)
# p1 <- p |> add_ground(var = "Pmask_dark_blue")
# p1 |> add_fill(var = "PRECAST_k16_nnSVG")
# 
# dev.off()
# 
# 

## 
pdf(here("plots", "VistoSeg", "nuclei density plots.pdf"), width = 10, height = 10)
for (sample_id in unique(colData(spe)$sample_id)){
  speb = spe[, which(spe$sample_id == sample_id)]
  den <- density(speb$Pmask_dark_blue)
  plot(den, frame = FALSE, col = "blue",main = sample_id)
}
dev.off()

colData(spe)$nuc_bin = 0
colData(spe)$nuc_bin[colData(spe)$Pmask_dark_blue>0.05] = 1
colData(spe)$nuc_bin[colData(spe)$Pmask_dark_blue>0.25] = 2
colData(spe)$nuc_bin = as.factor(colData(spe)$nuc_bin)
spe$nuc_bin = factor(spe$nuc_bin, levels = c("1" ,"2"))

pdf(here("plots", "VistoSeg", "escheR_3bin_emptyneuropilspots.pdf"), width = 10, height = 10)
for (sample_id in unique(colData(spe)$sample_id)){
speb = spe[, which(spe$sample_id == sample_id)]
p <- make_escheR(speb)
p1 <- p |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")
p1 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))
print(p1)
}
dev.off()

pdf(here("plots", "VistoSeg", "escheR_3bin_emptyneuropilspots.pdf"), width = 15, height = 15)

# for (sample_id in unique(spe$sample_id)){
# speb = spe[, which(spe$sample_id == sample_id)]; 
# p <- make_escheR(speb)
# p = p |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")+theme(legend.position = "none")
# p = p |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))
# print(p)
# }
# dev.off()

for (brain in unique(spe$brnum)){
  speb <- spe[, which(spe$brnum == brain)]
  samples <- unique(speb$sample_id)
  samples
  
   if (length(samples) == 2){
    spea = speb[, which(speb$sample_id == samples[1])]; p1 = make_escheR(spea)
    p1 = p1 |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")+theme(legend.position = "none")
    p1 = p1 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))
   
    spea = speb[, which(speb$sample_id == samples[2])]; p2 = make_escheR(spea)
    p2 = p2 |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")+ theme(legend.position = "none")
    p2 = p2 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))
    
    grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
    spea = speb[, which(speb$sample_id == samples[1])]; p1 = make_escheR(spea)
    p1 = p1 |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")+ theme(legend.position = "none")
    p1 = p1 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))
    
    spea = speb[, which(speb$sample_id == samples[2])]; p2 = make_escheR(spea)
    p2 = p2 |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")+ theme(legend.position = "none")
    p2 = p2 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))
    
    spea = speb[, which(speb$sample_id == samples[3])]; p3 = make_escheR(spea)
    p3 = p3 |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")+ theme(legend.position = "none")
    p3 = p3 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))
    
    grid.arrange(p1, p2, p3, nrow = 2)
    } else if (length(samples) == 4){
    spea = speb[, which(speb$sample_id == samples[1])]; p1 = make_escheR(spea)
    p1 = p1 |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")+ theme(legend.position = "none")
    p1 = p1|> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))
    
    spea = speb[, which(speb$sample_id == samples[2])]; p2 = make_escheR(spea)
    p2 = p2 |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")+ theme(legend.position = "none")
    p2 = p2|> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))
    
    spea = speb[, which(speb$sample_id == samples[3])]; p3 = make_escheR(spea)
    p3 = p3 |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")+ theme(legend.position = "none")
    p3 = p3|> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))

    spea = speb[, which(speb$sample_id == samples[4])]; p4 = make_escheR(spea)
    p4 = p4 |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")+ theme(legend.position = "none")
    p4 = p4|> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))
    
    grid.arrange(p1, p2, p3, p4, nrow = 2)
    } else if (length(samples) == 5){
    spea = speb[, which(speb$sample_id == samples[1])]; p1 = make_escheR(spea)
    p1 = p1 |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")+ theme(legend.position = "none")
    p1 = p1|> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))
    
    spea = speb[, which(speb$sample_id == samples[2])]; p2 = make_escheR(spea)
    p2 = p2 |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")+ theme(legend.position = "none")
    p2 = p2 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))
    
    spea = speb[, which(speb$sample_id == samples[3])]; p3 = make_escheR(spea)
    p3 = p3 |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")+ theme(legend.position = "none")
    p3 = p3|> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))
    
    spea = speb[, which(speb$sample_id == samples[4])]; p4 = make_escheR(spea)
    p4 = p4 |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")+ theme(legend.position = "none")
    p4 = p4 |> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))
    
    spea = speb[, which(speb$sample_id == samples[4])]; p5 = make_escheR(spea)
    p5 = p5 |> add_ground(var = "PRECAST_k18_nnSVG")+scale_colour_viridis_d(option = "H")+ theme(legend.position = "none")
    p5 = p5|> add_symbol(var = "nuc_bin", size = 0.4)+scale_shape_manual(values=c(20, 19))
    
    grid.arrange(p1, p2, p3, p4, p5, nrow = 2)}
}
dev.off()
